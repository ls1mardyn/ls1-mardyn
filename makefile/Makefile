# compiler configuration
CFG := gcc
#$(info $(CFG) set through $(origin CFG))

# machine type: SEQ or PAR
PARTYPE := SEQ
#$(info $(PARTYPE) set through $(origin PARTYPE))

# choose target: DEBUG|RELEASE
TARGET := DEBUG
#$(info $(TARGET) set through $(origin TARGET))


CURDIR := $(realpath $(PWD) )
ifeq ($(findstring src,$(CURDIR)),)
 $(warning ERROR: current directory $(CURDIR) seems not to be a source directory)
endif

SOURCES_COMMON = $(shell find ./ -name "*.cpp" | grep -v "parallel/")
SOURCES_SEQ = $(shell find parallel/ -name "*.cpp" | grep "Dummy")
SOURCES_PAR = $(shell find parallel/ -name "*.cpp" | grep -v "Dummy")
SOURCES = $(SOURCES_COMMON) $(SOURCES_$(PARTYPE))

OBJECTS_PAR = $(SOURCES_PAR:.cpp=.o)
OBJECTS = $(SOURCES:.cpp=.o)

# boost dependency
#BOOST_VER := boost_1_34_0

#LIBRARY_BOOST_LOCAL = $(shell find ../libs/$(BOOST_VER) -name 'libboost_program_options*.a' 2>/dev/null |grep release |grep threading-multi)
#LIBRARY_BOOST_SYSTEM = $(shell find /usr/lib/ -name 'libboost_program_options*.a' 2>/dev/null |grep -e '-mt' ; find /usr/local/lib/ -name 'libboost_program_options*.a' |grep -e '-mt')

# tinyXPATH dependency
#OBJECTS += $(patsubst %.cpp,%.o,$(shell find External/ticpp -name '*.cpp'))

INCLUDES = -I$(shell pwd)

#ifeq ($(strip $(LIBRARY_BOOST_LOCAL)),)
# $(info using BOOST library from system: $(LIBRARY_BOOST_SYSTEM))
# INCLUDES += -I$(shell dirname `echo $(LIBRARY_BOOST_SYSTEM) |sed -e 's/\ .*//g'`)
# OBJECTS += $(LIBRARY_BOOST_SYSTEM)
#else
# $(info using BOOST library from ls1 SVN: $(LIBRARY_BOOST_LOCAL))
# INCLUDES += -I../libs/$(BOOST_VER)
# OBJECTS += $(LIBRARY_BOOST_LOCAL)
#endif

#INCLUDES += -I$(PWD)/External/ticpp

BINARYbasename := MarDyn
BINARY := ${BINARYbasename}_$(shell svnversion)
BINARY := ${BINARY}_${PARTYPE}.${TARGET}

MAKEFILEDIR = $(dir $(realpath $(lastword $(MAKEFILE_LIST))))
CFGFILE := $(realpath $(MAKEFILEDIR)cfg/$(CFG).in )

# which compiler?
#include ../makefile/cfg/$(CFG).in
#ifndef (CFGFILE,)
# $(error ERROR: config file for $(CFG) not found - no include )
#else
# $(info including $(CFGFILE) )
 include $(CFGFILE)
#endif

CXX = $(CXX_$(PARTYPE))
CXXFLAGS = $(CXXFLAGS_$(PARTYPE)_$(TARGET)) -DMPICH_IGNORE_CXX_SEEK
LINKFLAGS = $(LINKFLAGS_$(PARTYPE)_$(TARGET)) -lrt

# only linker, domain decomposition and Simulation.cpp need MPI
$(BINARY): $(OBJECTS)
	@echo "building $(BINARY)"
	$(CXX) $(LINKFLAGS) -o $(BINARY) $(OBJECTS)
	$(shell if [ ! -f $(BINARYbasename) -o -L $(BINARYbasename) ]; then ln -sf $(BINARY) $(BINARYbasename); fi)

Simulation.o: Simulation.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(OBJECTS_PAR): %.o: %.cpp
	$(CXX) -c $(CXXFLAGS) $< -o $@

# all other files can always be compiled with g++
%.o: %.cpp
	$(CXX_SEQ) $(CXXFLAGS) -c $< -o $@

.PHONY: help all clean cleanall depend check boost cleanboost

help:
	@echo "configuration: $(CFG_INFO)"
	@echo "targets:"
	@echo "make all		build $(BINARY)"
	@echo "make clean	delete object files"
	@echo "make cleanall	delete object files and executable"
	@echo "make depend	generate makefile dependency list"
	@echo "make check	check for needed libraries"
	@echo "make boost	build bundled boost library"
	@echo "make cleanboost	clean bundled boost library"

all: check $(BINARY)

clean:
	find -name '*.o' -exec rm -r {} \;

cleanall: clean
	rm -f $(BINARY)
	if [ -L $(BINARYbasename) ]; then rm -f $(BINARYbasename); fi

depend:
	makedepend  -I$(INCLUDES) -- $(CXXFLAGS) -- $(SOURCES)

#check:
#ifeq ($(strip $(LIBRARY_BOOST_LOCAL)),)
#ifeq ($(strip $(LIBRARY_BOOST_SYSTEM)),)
#	@echo "ERROR: no boost library found. please run 'make boost' or install the boost 1.34.0 or newer package."
#else
#	@echo "SUCCESS: found libboost in $(shell dirname `echo $(LIBRARY_BOOST_SYSTEM) |sed -e 's/\ .*//g'`)"
#	@echo " includes: $(INCLUDES)"
#endif
#else
#	@echo "SUCCESS: found libboost in $(shell dirname `echo $(LIBRARY_BOOST_LOCAL) |sed -e 's/\ .*//g'`)"
#	@echo " includes: $(INCLUDES)"
#endif

boost:
	cd ../libs/$(BOOST_VER) && ./configure --with-libraries=program_options
	cd ../libs/$(BOOST_VER) && make

cleanboost:
	cd ../libs/$(BOOST_VER) && make distclean

# DO NOT DELETE THIS LINE -- make depend depends on it.
