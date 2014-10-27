# Set the task name
TASK = runasp

# Versions
VERSION = `python runasp.py --code-version`

# Uncomment the correct choice indicating either SKA or TST flight environment
FLIGHT_ENV = SKA

BIN = runasp
SHARE = runasp.py VERSION

include /proj/sot/ska/include/Makefile.FLIGHT

install:
#  Uncomment the lines which apply for this task
	mkdir -p $(INSTALL_BIN)
	mkdir -p $(INSTALL_SHARE)
#	mkdir -p $(INSTALL_DATA)
#	mkdir -p $(INSTALL_DOC)

	rsync --times --cvs-exclude $(BIN) $(INSTALL_BIN)/
	rsync --times --cvs-exclude $(SHARE) $(INSTALL_SHARE)/
#	rsync --times --cvs-exclude $(DATA) $(INSTALL_DATA)/
#	rsync --recursive --links --times -D --cvs-exclude $(DOC)/ $(INSTALL_DOC)/

