#  This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
#  It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
#  and by which you are legally bound while utilizing this file in any form or way.
#  There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#  You should have received a copy of that license along with this program;
#  if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
#  write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
#  Copyright (c) 2002-2003; all rights reserved unless otherwise stated.
# 
########### load platform dependent part
PLATFORM=$(shell uname -s)
include Makefile.$(PLATFORM)

########### list of source files
CSOURCES= main initialize global communication binary_file interaction_data \
	  verlet grid integrate cells ghosts forces rotation debug particle_data \
	  thermostat statistics statistics_chain energy pressure vmdsock imd \
	  p3m fft random blockfile blockfile_tcl polymer specfunc mmm1d tuning \
	  uwerr parser
CXXSOURCES=

LIBOBJECTS=c_blockfile.o

########### RULES
#################
TCLMD_CFLAGS=-DTCL_FILE_IO
OBJECTS=$(CSOURCES:%=%.o) $(CXXSOURCES:%=%.o)
CFILES=$(CSOURCES:=.c)
CXXFILES=$(CXXSOURCES:=.cc)

DOCFILES=$(shell ls doc/text/*.doc)

default: $(PLATFORM) $(PLATFORM)/Espresso $(PLATFORM)/libEspresso.a
all: $(PLATFORM) $(PLATFORM)/Espresso $(PLATFORM)/libEspresso.a

########### documentation
docu: doc/html/index.html

doc/html/index.html: $(DOCFILES) $(CFILES) $(CXXFILES)
	doxygen doxygen_config | grep -ve "^\(Generating\|Parsing\|Preprocessing\)"
#       (cd doc/latex; make)

########### output directory
$(PLATFORM):
	mkdir -p $(PLATFORM)

########### final target
$(PLATFORM)/Espresso: $(OBJECTS)
	(cd $(PLATFORM); $(LINK) $(LDFLAGS) -o Espresso $(OBJECTS) $(LDLIBS) $(STATIC_POSTLOAD))

$(PLATFORM)/libEspresso.a: $(LIBOBJECTS)
	(cd $(PLATFORM); ar -crs libEspresso.a $(LIBOBJECTS) )

########### clean
clean:
	rm -f *~
	(cd $(PLATFORM); rm -f $(OBJECTS) )
docclean:
	rm -rf doc/html/* doc/rtf/* doc/latex/* doc/man/*
mostclean: clean docclean
	rm -rf $(PLATFORM)

########### dependencies
dep: 
	$(MAKE) $(PLATFORM)
	rm -f $(PLATFORM)/.depend
	$(MAKE) $(PLATFORM)/.depend

$(PLATFORM)/.depend:
	mkdir -p $(PLATFORM)
	rm -f $@
	touch $@
	$(DEPEND) -f $@ -- $(CFLAGS) -- $(CFILES) $(CXXFILES) 2>/dev/null

include $(PLATFORM)/.depend

########## tests
test: $(PLATFORM)/Espresso
	cd testsuite; ./test.sh

########## implicit rules
vpath %.o  $(PLATFORM)

%.o: %.c
	$(CC) $(CFLAGS) $(TCLMD_CFLAGS) -c -o $(PLATFORM)/$@ $<

c_blockfile.o: blockfile.c
	$(CC) $(CFLAGS) -c -o $(PLATFORM)/$@ $<

%.o: %.cc
	$(CXX) $(CXXFLAGS) $(TCLMD_CFLAGS) -c -o $(PLATFORM)/$@ $<
