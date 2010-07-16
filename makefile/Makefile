# compiler configuration
CFG ?= gcc

# machine type: SEQ or PAR
PARTYPE ?= SEQ

# choose target: DEBUG|RELEASE
TARGET ?= RELEASE

# compile using depmode?
DEPS ?= 0

#build unit tests?
TESTS ?= 1

CURDIR ?= $(realpath $(PWD) )
ifeq ($(findstring src,$(CURDIR)),)
 $(warning ERROR: current directory $(CURDIR) seems not to be a source directory)
endif

#build name and include config file
MAKEFILEDIR := $(dir $(realpath $(lastword $(MAKEFILE_LIST))))
CFGFILE := $(realpath $(MAKEFILEDIR)cfg/$(CFG).in )

ifneq ($(CFGFILE),)
 $(info including $(CFGFILE) )
 include $(CFGFILE)
else
 $(error ERROR: config file for $(CFG) not found - no include )
endif

SOURCES_COMMON = $(shell find ./ -name "*.cpp" | grep -v "parallel/" | grep -v "/tests/")
SOURCES_SEQ = $(shell find parallel/ -name "*.cpp" | grep "Dummy" | grep -v "/tests/")
SOURCES_PAR = $(shell find parallel/ -name "*.cpp" | grep -v "Dummy" | grep -v "/tests/")
SOURCES = $(SOURCES_COMMON) $(SOURCES_$(PARTYPE)) 


OBJECTS = $(SOURCES:.cpp=.o)
# boost dependency
#BOOST_VER := boost_1_34_0

#LIBRARY_BOOST_LOCAL = $(shell find ../libs/$(BOOST_VER) -name 'libboost_program_options*.a' 2>/dev/null |grep release |grep threading-multi)
#LIBRARY_BOOST_SYSTEM = $(shell find /usr/lib/ -name 'libboost_program_options*.a' 2>/dev/null |grep -e '-mt' ; find /usr/local/lib/ -name 'libboost_program_options*.a' |grep -e '-mt')

# tinyXPATH dependency
#OBJECTS += $(patsubst %.cpp,%.o,$(shell find External/ticpp -name '*.cpp'))

#INCLUDES = -I$(shell pwd)
INCLUDES = -I$(CURDIR) -I$(CURDIR)/External $(CPPUNIT_INCLUDES)

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

BINARYbasename = MarDyn
ifneq ($(realpath .svn),)
  $(info Building from a svn version.)
  svnVERSION = $(shell svnversion | sed -e "s/:/-/")
  ifneq ($(svnVERSION),)
    BINARY = $(BINARYbasename)_$(svnVERSION).$(PARTYPE)_$(TARGET)
  else
    BINARY = $(BINARYbasename).$(PARTYPE)_$(TARGET)
  endif
else
  BINARY = $(BINARYbasename).$(PARTYPE)_$(TARGET)
endif

CXX = $(CXX_$(PARTYPE))
CXXFLAGS = $(CXXFLAGS_$(PARTYPE)_$(TARGET)) -DMPICH_IGNORE_CXX_SEEK
LINKFLAGS = $(LINKFLAGS_$(PARTYPE)_$(TARGET)) $(CPPUNIT_LINK_ARGS) -lrt

ifeq ($(DEPS),1)
DEPCOMP = depmode=$(DEPMODE) source=$< object=$@ $(MAKEFILEDIR)depcomp
else
DEPCOMP =
endif

$(BINARY): .deps $(OBJECTS) cppunit
	@echo "building $(BINARY)"
	$(CXX) $(LINKFLAGS) -o $@ $(OBJECTS)
	$(shell if [ ! -f $(BINARYbasename) -o -L $(BINARYbasename) ]; then ln -sf $(BINARY) $(BINARYbasename); fi)

%.o: %.cpp
	$(DEPCOMP) $(CXX) $(CXXFLAGS) -c $< -o $@

.PHONY: help all clean cleanall depend check boost cleanboost

help:
	@echo "configuration: $(CFG_INFO)"
	@echo
	@echo "options:"
	@echo "make CFG=...                 choose compiler (see make cfg_list)"
	@echo "make PARTYPE=SEQ | PAR       sequential or parallel (MPI) build"
	@echo "make TARGET=DEBUG | RELEASE  debugging / non optimized build"
	@echo "make DEPS=0 | 1              compile using dependancy tracking"
	@echo "make TESTS=0 | 1     compile with or without unit tests (and build cppunit if neccessary)"
	@echo
	@echo "targets:"
	@echo "make all         build $(BINARY)"
	@echo "make cfg_list    list available configurations"
	@echo "make clean       delete object files"
	@echo "make cleanall    delete object files and executable"
	@echo "make depend      generate makefile dependency list"
	#@echo "make check       check for needed libraries"
	#@echo "make boost       build bundled boost library"
	#@echo "make cleanboost  clean bundled boost library"

all: check $(BINARY)

test: 
	./$(BINARY) -t

#list available configurations
cfg_list:
	@ls -1 ../makefile/cfg | sed -e 's/.in//'

clean:
	find -type f -name '*.o' -delete
	find -path '*/.deps/*' -type f -delete
	find -type d -name .deps -delete

clear: cleanall

cleanall: clean cppunit_clean
	rm -f $(BINARY)
	find -type f -name '*~' -delete
	if [ -L $(BINARYbasename) ]; then rm -f $(BINARYbasename); fi

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

# create .deps subdirectories containing dependency information
.deps:
ifeq ($(DEPS),1)
	$(shell find -type d ! -name .deps -execdir mkdir -p {}/.deps \;)
endif

# include potentially generated dependency information
-include $(shell find -path '*/.deps/*' -type f)
  
#include makefile for cppunit tests
CPPUNIT_MAKEFILE := $(realpath $(MAKEFILEDIR)cppunit.mk )
include $(CPPUNIT_MAKEFILE)

