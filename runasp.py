#!/usr/bin/env python
import os
import sys
import re
from glob import glob
import csv
import gzip
import shutil
from itertools import izip, count
import numpy as np

import Ska.arc5gl
from Ska.Shell import getenv, bash, tcsh_shell, ShellError
import pyyaks.logger
from astropy.io import fits

_versionfile = os.path.join(os.path.dirname(__file__), 'VERSION')
VERSION = open(_versionfile).read().strip()


def get_options():
    from optparse import OptionParser
    parser = OptionParser()
    parser.set_defaults()
    parser.add_option("--obsid",
                      type='int',
                      help="obsid to process")
    parser.add_option("--version",
                      help="version of products")
    parser.add_option("--revision",
                      default=1,
                      type="int",
                      help="integer revision label for output")
    parser.add_option("--label",
                      help="output processing dir toplevel name")
    parser.add_option("--skip-slot",
                      type='int',
                      action='append',
                      help="slots to skip in processing")
    parser.add_option('--skip-slot-method',
                      default='telem',
                      help="method to cut slots")
    parser.add_option("--range",
                      help="processing range specifier")
    parser.add_option("--dir",
                      default="pipeline_out",
                      help="directory for telem fetching")
    parser.add_option("--pipe-stop-before",
                      help="stop pipeline before specified pipe")
    parser.add_option("--pipe-start-at",
                      help="start pipeline at specified pipe")
    parser.add_option("--code-version",
                      action='store_true',
                      help="return version of the runasp tool")
    parser.add_option("--log-level",
                      type='int',
                      default=20,
                      help="Log level, default=20 (10=debug, 15=verbose, 20=info, 100=quiet)")
    parser.add_option('--log-file',
                      default='pipe.log',
                      help='Log file (default=pipe.log)')
    opt, args = parser.parse_args()
    return opt, args


pcad0_files = ['pcad0/pcad*_%d_*fits*' % x for x in [3, 5, 7, 8, 14, 15]]

pipe_config = dict(pipe_ped='asp_l1_std',
                   infiles=pcad0_files
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
cai_override = {'obs_id': 'i',
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


def dir_setup(dir, istart, label=None, inplace=False, rev=1):
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
    if len(indirs) and not inplace:
        rev = len(indirs) + 1
        indir = os.path.join(workdir, "in%d" % rev)
        if os.path.exists(indir):
            raise ValueError("Bad in directory sequence (%s exists)"
                             % indir)
    if not os.path.exists(indir):
        logger.info('Making input directory {}'.format(indir))
        os.makedirs(indir)
    outdir = os.path.join(workdir, "out%d" % rev)
    if os.path.exists(outdir) and not inplace:
        raise ValueError("Bad in directory sequence (%s exists)"
                         % outdir)
    if not os.path.exists(outdir):
        os.makedirs(outdir)
        logger.info('Making output directory {}'.format(outdir))
    return workdir, indir, outdir


def link_files(dir, indir, outdir, istart, istop, obiroot, skip_slot=None):
    """
    Creates symbolic links from the specified indir to the available telemetry.
    Fits files are only linked in if their header time keywords are relevant.
    ACA0 image files may be skipped if the slot is in skip_slot list.
    obspars must be for the correct obi.
    """
    dirmap = dict(infiles=indir,
                  outfiles=outdir)
    for filetype in ['infiles', 'outfiles']:
        ldir = dirmap[filetype]
        for fileglob in pipe_config[filetype]:
            match = glob(os.path.join(opt.dir, fileglob))
            for mfile in match:
                fitsmatch = re.match('.*fits', mfile)
                if fitsmatch:
                    header = fits.getheader(mfile)
                    if ((istart >= header['tstop'])
                            or (istop <= header['tstart'])):
                        logger.verbose("skipping file out of timerange {}".format(mfile))
                        continue
                    aca0 = re.search('aca.*_(\d)_img0', mfile)
                    if skip_slot and aca0:
                        aca_file_slot = int(aca0.group(1))
                        if aca_file_slot in skip_slot:
                            logger.verbose("skipping slot file on {}".format(mfile))
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


def get_range_ai(ai_cmds, proc_range):
    """
    Limit the array of "todo" aspect to those requested in the --range
    specifier
    """
    # if no --range, return all
    if not proc_range:
        return ai_cmds
    # if a single integer, return on that aspect interval
    intmatch = re.match('^(\d+)$', proc_range)
    if intmatch:
        interv = int(intmatch.group(1))
        return [ai_cmds[int(intmatch.group(1))]]
    # if of the form 0:1, return that range of intervals
    # (python form, not inclusive)
    imatch = re.match('^(\d+):(\d+)$', proc_range)
    if imatch:
        return ai_cmds[int(imatch.group(1)):int(imatch.group(2))]
    # if of the form 1: , return range 1 -> end
    omatch = re.match('^(\d+):$', proc_range)
    if omatch:
        return ai_cmds[int(omatch.group(1)):]
    # if of the form 0:+3000, find a tstop corresponding
    # to tstart of aspect interval 0 plus 3000 seconds
    tmatch = re.match('^(\d+):\+(\d+)$', proc_range)
    if tmatch:
        # get n seconds of specified interval
        interv = int(tmatch.group(1))
        seconds = int(tmatch.group(2))
        tstart = ai_cmds[interv]['istart']
        tstop = tstart + seconds
        cut_ai = ai_cmds[interv].copy()
        cut_ai['istop'] = tstop
        print "attempted to process first %d sec of ai %d" % (
            seconds, interv)
        return [cut_ai]


def cut_stars(ai):
    starfiles = glob(os.path.join(ai['outdir'],
                                  "*stars.txt"))
    shutil.copy(starfiles[0], starfiles[0] + ".orig")
    starlines = open(starfiles[0]).read().split("\n")
    for slot in ai['skip_slot']:
        starlines = [i for i in starlines
                     if not re.match("^\s+{}\s+1.*".format(slot), i)]
    logger.info('Cutting stars by updating {}'.format(starfiles[0]))
    with open(starfiles[0], "w") as newlist:
        newlist.write("\n".join(starlines))


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


def run_ai(ais):
    """
    Run aspect pipeline 'flt_run_pipe' over the aspect intervals described
    in the list of dictionaries passed as an argument
    """
    ascds_env = getenv('source /home/ascds/.ascrc -r release', shell='tcsh')
    ocat_env = getenv(
        'source /proj/sot/ska/data/aspect_authorization/set_ascds_ocat_vars.csh',
        shell='tcsh')
    for var in ['ASCDS_OCAT_UNAME', 'ASCDS_OCAT_SERVER', 'ASCDS_OCAT_PWORD']:
        ascds_env[var] = ocat_env[var]

    logger_fh = FilelikeLogger(logger)

    for ai in ais:
        pipe_cmd = 'flt_run_pipe -r {root} -i {indir} -o {outdir} \
-t {pipe_ped} \
-a "INTERVAL_START"={istart} \
-a "INTERVAL_STOP"={istop} \
-a obiroot={obiroot} \
-a revision=1 '.format(**ai)
        if 'pipe_start_at' in ai:
            pipe_cmd = pipe_cmd + " -s {}".format(ai['pipe_start_at'])
        if 'pipe_stop_before' in ai:
            pipe_cmd = pipe_cmd + " -S {}".format(ai['pipe_stop_before'])
        if 'skip_slot' in ai:
            try:
                tcsh_shell(pipe_cmd + " -S check_star_data",
                           env=ascds_env,
                           logfile=logger_fh)
            except ShellError as sherr:
                # if shell error, just check to see if get_star_data completed successfully
                loglines = open(logger_fh.filename).read()
                if not re.search("get_star_data completed successfully", loglines):
                    raise ShellError(sherr)
            cut_stars(ai)
            tcsh_shell(pipe_cmd + " -s check_star_data",
                       env=ascds_env,
                       logfile=logger_fh)
        else:
            logger.info('Running pipe command {}'.format(pipe_cmd))
            tcsh_shell(pipe_cmd,
                       env=ascds_env,
                       logfile=logger_fh)


def mock_cai_file(opt):
    """
    Mock up a Constant Aspect Interval file
    based on the Kalman intervals in the aiprops files
    """
    aiprops_files = glob(os.path.join(opt.dir, "asp05/*aipr*"))
    obs_records = []
    colnames = ['obsid', 'pcad_mode', 'aspect_mode',
                'start_time', 'stop_time']
    for ai_file in aiprops_files:
        tbl = fits.open(ai_file)[1].data
        for row in tbl[tbl['obsid'] == opt.obsid]:
            obs_records.append([row[s] for s in colnames])

    ai_rec = np.rec.fromrecords(obs_records, names=colnames)
    acq = np.flatnonzero(ai_rec['aspect_mode'] == 'ACQUISITION')
    gui = np.flatnonzero(ai_rec['aspect_mode'] == 'GUIDE')
    if len(acq) and len(gui):
        kalman_mask = ((np.arange(len(ai_rec)) > acq)
                       & (np.arange(len(ai_rec)) > gui)
                       & (ai_rec['pcad_mode'] == 'NPNT')
                       & (ai_rec['aspect_mode'] == 'KALMAN'))
        kalman = ai_rec[kalman_mask]
    if not len(kalman):
        raise ValueError("no kalman intervals in aiprops")
    kalman_intervals = [dict((col, kalman[0][col]) for col in colnames)]

    if len(kalman) > 1:
        for k, idx in izip(kalman[1:], count(1)):
            if k['start_time'] == kalman[idx - 1]['stop_time']:
                kalman_intervals[-1].update(dict(stop_time=k['stop_time']))
            else:
                kalman_intervals.append(k)
    obspar_files = glob(os.path.join(opt.dir, "obspar/*.par*"))
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
    cai_file = os.path.join(opt.dir, "asp05", "pcadf%05d_000N%03d_cai0a.par"
                            % (int(obspar['obs_id']), obspar['revision']))
    fh = open(cai_file, 'w')
    fh.write(cai_text)
    fh.close()


def main(opt):
    # get files
    if opt.obsid:
        logger.info('Opening connection to archive server')
        arc5 = Ska.arc5gl.Arc5gl()
        for (prod, query) in pipe_config['archfiles']:
            proddir = os.path.join(opt.dir, prod)
            if not os.path.exists(proddir):
                logger.info('Creating directory {}'.format(proddir))
                os.makedirs(proddir)
            else:
                logger.info('Skipping directory {}: exists'.format(proddir))
                continue
            obsid = opt.obsid
            if logger.level < 20:
                arc5.echo = True
            arc5.sendline("cd %s" % os.path.abspath(proddir))
            arc5.sendline("obsid=%d" % int(obsid))
            if opt.version is not None:
                arc5.sendline("version=%s" % opt.version)
            logger.info('Sending "get {}"'.format(query))
            arc5.sendline("get %s" % query)
            gotfiles = glob(os.path.join(proddir, "*"))
            if not len(gotfiles):
                os.rmdir(proddir)
        del arc5

    caiprops_files = glob(os.path.join(opt.dir, "asp05",
                                       "pcad*cai0a.par*"))
    if not len(caiprops_files):
        if os.path.exists(os.path.join(opt.dir, "asp05")):
            mock_cai_file(opt)
        else:
            raise ValueError

    # check files
    for filetype in ['infiles', 'outfiles']:
        for fileglob in pipe_config[filetype]:
            match = glob(os.path.join(opt.dir, fileglob))
            if not len(match):
                raise ValueError("No files found for glob %s"
                                 % fileglob)
            for mfile in match:
                if re.match(".*\.gz", mfile):
                    logger.verbose('Unzipping {}'.format(mfile))
                    bash("gunzip -f %s" % os.path.abspath(mfile))

    # reset this to get unzipped names
    caiprops_files = glob(os.path.join(opt.dir, "asp05",
                                       "pcad*cai0a.par*"))

    # constant aspect interval files
    obi = {}
    if len(caiprops_files):
        for cai_file in caiprops_files:
            cai = get_par(cai_file, cai_override)
            if not cai['obi_num'] in obi:
                obi[cai['obi_num']] = {}
            interval = 0
            while ("istart_%d" % interval in cai
                   and "istop_%d" % interval in cai):
                obi[cai['obi_num']][interval] = \
                    {'istart': cai['istart_%d' % interval],
                     'istop':  cai['istop_%d' % interval]}
                interval += 1

    ai_cmds = []
    for obi_num in obi:
        # read possible obspars
        obspar_files = glob(os.path.join(opt.dir, "obspar/*.par"))
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
            # Work inplace?  Only if told to start at a specific pipe
            inplace = False
            if opt.pipe_start_at:
                inplace = True
            # directory setup
            workdir, indir, outdir = dir_setup(opt.dir,
                                               int(istart),
                                               label=opt.label,
                                               inplace=inplace,
                                               rev=opt.revision)

            # if skipping the slot by chucking the telem
            telem_skip_slot = []
            process_skip_slot = []
            if opt.skip_slot_method == 'telem':
                telem_skip_slot = opt.skip_slot
            else:
                process_skip_slot = opt.skip_slot

            # link relevant files
            link_files(opt.dir, indir, outdir, istart, istop,
                       obiroot, telem_skip_slot)

            # make list files
            make_list_files(opt.dir, indir, outdir, root)

            # spec
            cmd = dict(dir=os.path.abspath(opt.dir),
                       obi=obi_num,
                       indir=indir,
                       outdir=outdir,
                       root=root,
                       pipe_ped="%s.ped" % pipe_config['pipe_ped'],
                       istart=istart,
                       istop=istop,
                       obiroot=obiroot,
                       log="%s_f%09d.log" % (pipe_config['pipe_ped'], istart))
            if len(process_skip_slot):
                cmd['skip_slot'] = process_skip_slot
            if opt.pipe_start_at:
                cmd['pipe_start_at'] = opt.pipe_start_at
            if opt.pipe_stop_before:
                cmd['pipe_stop_before'] = opt.pipe_stop_before

            ai_cmds.append(cmd)

    range_ais = get_range_ai(ai_cmds, opt.range)
    run_ai(range_ais)

if __name__ == '__main__':
    global logger
    opt, args = get_options()
    if opt.code_version:
        print VERSION
        sys.exit(0)
    if opt.pipe_stop_before == 'create_props_files':
        print """There is a known bug in flt_pctr that can prevent
a stop at 'create_props_files'.  Choose a different pipe stop."""
        sys.exit(0)
    logger = pyyaks.logger.get_logger(name='runasp', level=opt.log_level, filename=opt.log_file,
                                      filelevel=15, format="%(asctime)s %(message)s")
    main(opt)
