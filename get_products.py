#!/proj/sot/ska/bin/python

import os
import re
from glob import glob
import csv
import gzip
import numpy as np
import json

import Ska.arc5gl
from Ska.Shell import bash
import pyyaks.logger
from astropy.io import fits
from astropy.table import Table


def get_options():
    from optparse import OptionParser
    parser = OptionParser()
    parser.set_defaults()
    parser.add_option("--obsid",
                      type='int',
                      help="obsid to process")
    parser.add_option("--outfile",
                      default="interval_spec.json")
    parser.add_option("--version")
    parser.add_option("--label",
                      help="output processing dir toplevel name")
    parser.add_option("--dir",
                      default="pipeline_out",
                      help="directory for telem fetching")
    parser.add_option("--log-level",
                      type='int',
                      default=20,
                      help="Log level, default=20 (10=debug, 15=verbose, 20=info, 100=quiet)")
    parser.add_option('--log-file',
                      default='get_products.log',
                      help='Log file (default=pipe.log)')
    opt, args = parser.parse_args()
    return opt, args


PCAD0_FILES = ['pcad0/pcad*_%d_*fits*' % x for x in [3, 5, 7, 8, 14, 15]]

PIPE_CONFIG = dict(pipe_ped='asp_l1_std',
                   infiles=PCAD0_FILES
                   + ['asp05/pcad*aipr0a.fits*',
                      'asp05/pcad*cai0a.*',
                      'aca0/*',
                      'sim05/sim*coor0a*',
                      'obspar/axaf*obs0a.par*',
                      'acis2eng/acis*fits*',
                      'obc4eng/obc*fits*'],
                   outfiles=[],
                   archfiles=[('aca0', 'aca0'),
                              ('pcad0', 'pcad0'),
                              ('sim05', 'sim05'),
                              ('asp05', 'asp05'),
                              ('obspar', 'obspar'),
                              ('acis2eng', 'acis_eng_0{acis2eng}'),
                              ('obc4eng', 'obc_eng_0{obc4eng}')])

# the cai files have non-compliant IRAF par files.
# override the data types
CAI_OVERRIDE = {'obs_id': 'i',
                'obi_num': 'i',
                'ascdsver': 's',
                'num_cai': 'i'}


def parse_obspar(file, override=None):
# borrowed from telem_archive ... override is new
    convert = {'i': int,
               'r': float,
               's': str}
    try:
        lines = gzip.open(file).readlines()
    except IOError:
        lines = open(file).readlines()
    obs_read = csv.DictReader(lines,
                              fieldnames=('name', 'type', 'hidden', 'value',
                                          'def1', 'def2', 'descr'),
                              dialect='excel')

    for row in obs_read:
        if override and (row['name'] in override):
            row['value'] = convert[override[row['name']]](row['value'])
        else:
            row['value'] = convert[row['type']](row['value'])
        row['name'] = row['name'].replace('-', '_')
        yield row

    return


def get_par(parfile, override=None):
    """
    Read an IRAF-style par file and return a dictionary.
    """
    par = {}
    for row in parse_obspar(parfile, override):
        par.update({row['name']: row['value']})
    return par


def get_obspar(obsparfile):
    """Get the obspar for obsid starting at tstart.  Return as a dict."""
    obspar = {'num_ccd_on': 0}
    for row in parse_obspar(obsparfile):
        obspar.update({row['name']: row['value']})
        if re.match(r'^ccd[is]\d_on$', row['name']) and row['value'] == 'Y':
            obspar['num_ccd_on'] += 1
    return obspar


def dir_setup(dir, istart, label=None, rev=1):
    """
    Makes

      - a directory to contain processing for each aspect interval (AI)
      - at least one input directory in each AI directory which will be
        populated with telemetry or links to telemetry for processing.
      - at least one output directory in each AI directory

    The input and output directories have trailing integer revision numbers
    that are incremented in successive processing attempts.
    """
    if label is None:
        workdir = os.path.join(dir,
                               'ASP_L1_STD_%09d' % istart)
    else:
        workdir = os.path.join(dir, label)
    if not os.path.exists(workdir):
        logger.info('Making working directory {}'.format(workdir))
        os.makedirs(workdir)
    else:
        logger.info('Using existing working directory {}'.format(workdir))

    indir = os.path.join(workdir, "in%d" % rev)
    indirs = glob(os.path.join(workdir, "in*"))
    if not os.path.exists(indir):
        logger.info('Making input directory {}'.format(indir))
        os.makedirs(indir)
    outdir = os.path.join(workdir, "out%d" % rev)
    if not os.path.exists(outdir):
        os.makedirs(outdir)
        logger.info('Making output directory {}'.format(outdir))
    return workdir, indir, outdir


def link_files(outdir, indir, rundir, istart, istop, obiroot):
    """
    Creates symbolic links from the specified indir to the available telemetry.
    Fits files are only linked in if their header time keywords are relevant.

    obspars must be for the correct obi.
    """
    dirmap = dict(infiles=indir,
                  outfiles=rundir)
    for filetype in ['infiles', 'outfiles']:
        ldir = dirmap[filetype]
        for fileglob in PIPE_CONFIG[filetype]:
            match = glob(os.path.join(outdir, fileglob))
            for mfile in match:
                fitsmatch = re.match('.*fits', mfile)
                if fitsmatch:
                    header = fits.getheader(mfile)
                    if ((istart >= header['tstop'])
                            or (istop <= header['tstart'])):
                        logger.verbose("skipping file out of timerange {}".format(mfile))
                        continue
                obsparmatch = re.match('.*obs0a\.par(\.gz)?', mfile)
                if obsparmatch:
                    obimatch = re.match('.*axaf%s_obs0a\.par(\.gz)?' % obiroot, mfile)
                    if not obimatch:
                        logger.verbose("skipping obspar for different obi")
                        continue
                if not os.path.exists(os.path.join(ldir, os.path.basename(mfile))):
                    logger.info("ln -s {} {}".format(os.path.relpath(mfile, ldir), ldir))
                    bash("ln -s %s %s" % (os.path.relpath(mfile, ldir), ldir))


def make_list_files(dir, indir, outdir, root):
    """
    Create .lis files for the pipeline.
    """
    # remove any present already
    inlists = glob(os.path.join(indir, "*.lis"))
    [os.unlink(x) for x in inlists]
    # remove any present already
    outlists = glob(os.path.join(outdir, "*.lis"))
    [os.unlink(x) for x in outlists]

    for listend, listglob in (('sim.lis', 'sim*coor0a*fits*'),
                              ('pcad.lis', 'pcad*eng0*fits*'),
                              ('acis.lis', 'acis*eng0*fits*'),
                              ('obc.lis', 'obc*eng0*fits*')):
        filename = os.path.join(indir, "{root}_{listend}".format(root=root, listend=listend))
        logger.info('Writing list file {}'.format(filename))
        with open(filename, 'w') as lfile:
            sglob = sorted(glob(os.path.join(indir, listglob)))
            lfile.write("\n".join([os.path.basename(x) for x in sglob]))

    # aca0
    filename = os.path.join(indir, "{root}_tel.lis".format(root=root))
    logger.info('Writing list file {}'.format(filename))
    with open(filename, 'w') as lfile:
        for slot in [3, 4, 5, 6, 7, 0, 1, 2]:
            sglob = sorted(glob(os.path.join(indir, 'aca*_%d_*0.fits*' % slot)))
            telem_lines = '\n'.join([os.path.basename(x) for x in sglob])
            lfile.write(telem_lines)
            lfile.write("\n")

    # pcad adat in outdir if present
    filename = os.path.join(outdir, "{root}_dat.lis".format(root=root))
    logger.info('Writing list file {}'.format(filename))
    with open(filename, 'w') as lfile:
        sglob = sorted(glob(os.path.join(indir, 'pcad*adat*fits*')))
        lfile.write('\n'.join([os.path.basename(x) for x in sglob]))



class FilelikeLogger(object):
    """
    Make logger object look a bit file-like for writing
    """
    def __init__(self, logger):
        self.logger = logger
        for fh in self.logger.handlers:
            try:
                self.filename = fh.baseFilename
                break
            except AttributeError:
                pass

    def write(self, value):
        self.logger.info(value.rstrip())

    def flush(self):
        for fh in self.logger.handlers:
            fh.flush()

    def close(self):
        pass



def mock_cai_file(obsid, outdir):
    """
    Mock up a Constant Aspect Interval file
    based on the Kalman intervals in the aiprops files
    """
    aiprops_files = glob(os.path.join(outdir, "asp05/*aipr*"))
    obs_records = []
    colnames = ['obsid', 'pcad_mode', 'aspect_mode',
                'start_time', 'stop_time']
    for ai_file in aiprops_files:
        tbl = fits.open(ai_file)[1].data
        for row in tbl[tbl['obsid'] == obsid]:
            obs_records.append([row[s] for s in colnames])

    ai_rec = np.rec.fromrecords(obs_records, names=colnames)
    acq = np.flatnonzero(ai_rec['aspect_mode'] == 'ACQUISITION')
    gui = np.flatnonzero(ai_rec['aspect_mode'] == 'GUIDE')
    if len(acq) and len(gui):
        kalman_mask = ((np.arange(len(ai_rec)) > np.max(acq))
                       & (np.arange(len(ai_rec)) > np.max(gui))
                       & (ai_rec['pcad_mode'] == 'NPNT')
                       & (ai_rec['aspect_mode'] == 'KALMAN'))
        kalman = Table(ai_rec[kalman_mask])
    if not len(kalman):
        raise ValueError("no kalman intervals in aiprops")
    kalman_intervals = [{col: kalman[0][col] for col in colnames}]
    # Compress any other intervals if last saved stop matches current start
    if len(kalman) > 1:
        for idx in range(1, len(kalman)):
            if np.abs(kalman[idx]['start_time'] - kalman_intervals[-1]['stop_time']) < .1:
                kalman_intervals[-1]['stop_time'] = kalman[idx]['stop_time']
            else:
                kalman_intervals.append({col: kalman[idx][col] for col in colnames})
    obspar_files = glob(os.path.join(outdir, "obspar/*.par*"))
    # just quit if there are multiple obspars
    if len(obspar_files) > 1:
        raise ValueError("Multiple obspars; unsure how to proceed.")
    obspar = get_obspar(obspar_files[0])
#    cai_text = "blah"
    cai_text = ('''obs_id,r,h,%(obs_id)s,,,""
obi_num,r,h,%(obi_num)d,,,""
start_time,r,h,%(tstart)f,,,""
stop_time,r,h,%(tstop)f,,,""
ascdsver,r,h,%(ascdsver)s,,,""
''' % obspar)
    cai_text = cai_text + 'num_cai,r,h,{},,,""\n'.format(len(kalman))
    for kal_int in kalman_intervals:
        cai_text = cai_text + ('''
istart_0,r,h,{},,,""
istop_0,r,h,{},,,""
'''.format(kal_int['start_time'], kal_int['stop_time']))
    cai_file = os.path.join(outdir, "asp05", "pcadf%05d_000N%03d_cai0a.par"
                            % (int(obspar['obs_id']), obspar['revision']))
    fh = open(cai_file, 'w')
    fh.write(cai_text)
    fh.close()


def get_archive_files(obsid, outdir, version=None):
    logger.info('Opening connection to archive server')
    arc5 = Ska.arc5gl.Arc5gl()
    for (prod, query) in PIPE_CONFIG['archfiles']:
        proddir = os.path.join(outdir, prod)
        if not os.path.exists(proddir):
            os.makedirs(proddir)
        if logger.level < 20:
            arc5.echo = True
        arc5.sendline("cd %s" % os.path.abspath(proddir))
        arc5.sendline("obsid=%d" % int(obsid))
        if version is not None:
            arc5.sendline("version=%s" % version)
        logger.info('Sending "get {}"'.format(query))
        arc5.sendline("get %s" % query)
        gotfiles = glob(os.path.join(proddir, "*"))
        if not len(gotfiles):
            os.rmdir(proddir)
    del arc5


def get_products(obsid, outdir, label, version):

    get_archive_files(obsid, outdir, version=version)
    caiprops_files = glob(os.path.join(outdir, "asp05",
                                       "pcad*{}*cai0a.par*".format(int(obsid))))
    if not len(caiprops_files):
        if os.path.exists(os.path.join(outdir, "asp05")):
            mock_cai_file(obsid, outdir)
        else:
            raise ValueError

    # check files
    for filetype in ['infiles', 'outfiles']:
        for fileglob in PIPE_CONFIG[filetype]:
            match = glob(os.path.join(outdir, fileglob))
            if not len(match):
                raise ValueError("No files found for glob %s"
                                 % fileglob)
            for mfile in match:
                if re.match(".*\.gz", mfile):
                    logger.verbose('Unzipping {}'.format(mfile))
                    bash("gunzip -f %s" % os.path.abspath(mfile))
    # reset this to get unzipped names
    caiprops_files = glob(os.path.join(outdir, "asp05",
                                       "pcad*{}*cai0a.par*".format(int(obsid))))

    # Read the constant aspect interval files to define the aspect intervals
    obi = {}
    if len(caiprops_files):
        for cai_file in caiprops_files:
            cai = get_par(cai_file, CAI_OVERRIDE)
            if not cai['obi_num'] in obi:
                obi[cai['obi_num']] = {}
            interval = 0
            while ("istart_%d" % interval in cai
                   and "istop_%d" % interval in cai):
                obi[cai['obi_num']][interval] = \
                    {'istart': cai['istart_%d' % interval],
                     'istop':  cai['istop_%d' % interval]}
                interval += 1

    # Set up aspect intervals (should really just be one for most observations)
    ai_cmds = []
    for obi_num in obi:
        # Read possible obspars
        obspar_files = glob(os.path.join(outdir, "obspar/*.par"))
        for ofile in obspar_files:
            obspar = get_obspar(ofile)
            if obspar['obi_num'] == obi_num:
                obsmatch = re.search('axaf(.+)_obs0a\.par', ofile)
                obiroot = obsmatch.group(1)
        if not obiroot:
            raise ValueError("no obspar for obi %d" % obi_num)

        for ai_num in obi[obi_num]:
            aspect_interval = obi[obi_num][ai_num]
            istart = aspect_interval['istart']
            istop = aspect_interval['istop']
            root = "f%09d" % istart
            # directory setup
            workdir, indir, rundir = dir_setup(outdir,
                                               int(istart),
                                               label=label)

            # link relevant files
            link_files(outdir, indir, rundir, istart, istop, obiroot)

            # make list files
            make_list_files(outdir, indir, rundir, root)

            # spec
            cmd = dict(outdir=os.path.abspath(outdir),
                       obi=obi_num,
                       indir=indir,
                       rundir=rundir,
                       root=root,
                       pipe_ped="%s.ped" % PIPE_CONFIG['pipe_ped'],
                       istart=istart,
                       istop=istop,
                       obiroot=obiroot,
                       log="%s_f%09d.log" % (PIPE_CONFIG['pipe_ped'], istart))
            ai_cmds.append(cmd)
    return ai_cmds

if __name__ == '__main__':
    global logger
    opt, args = get_options()
    logger = pyyaks.logger.get_logger(name='get_products.py', level=opt.log_level, filename=opt.log_file,
                                      filelevel=15, format="%(asctime)s %(message)s")
    ais = get_products(obsid=opt.obsid, outdir=opt.dir, label=opt.label, version=opt.version)
    json.dump(ais, open(opt.outfile, "w"), indent=4, sort_keys=True)
    raise ValueError
