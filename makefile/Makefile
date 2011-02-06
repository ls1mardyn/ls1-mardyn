# compiler configuration
# show list of all available configurations with make cfg_list
CFG ?= gcc

# machine type: SEQ or PAR
PARTYPE ?= SEQ

# choose target: DEBUG|RELEASE
# some configurations have additional targets, see the respective configuration files for details
TARGET ?= RELEASE

# compile using depmode?
DEPS ?= 0

#build unit tests?
TESTS ?= 0

VTK ?= 0

LCOV?= 0

CURDIR ?= $(realpath $(PWD) )
ifeq ($(findstring src,$(CURDIR)),)
  $(error ERROR: current directory $(CURDIR) seems not to be a source directory)
endif

#build name and include config file
MAKEFILEDIR := $(dir $(realpath $(lastword $(MAKEFILE_LIST))))
CFGFILE := $(realpath $(MAKEFILEDIR)cfg/$(CFG).in )

ifneq ($(CFGFILE),)
  $(info # Selected configuration: $(CFG) )
  $(info # Including file $(CFGFILE) )
  include $(CFGFILE)
else
  $(error ERROR: config file for $(CFG) not found - no include )
endif

SOURCES_COMMON = $(shell find ./ -name "*.cpp" | grep -v "parallel/" | grep -v "/tests/" | grep -v "/vtk/")
SOURCES_SEQ = $(shell find parallel/ -name "*.cpp" | grep "Dummy" | grep -v "/tests/" | grep -v "/vtk/")
SOURCES_PAR = $(shell find parallel/ -name "*.cpp" | grep -v "Dummy" | grep -v "/tests/" | grep -v "/vtk/")
SOURCES = $(SOURCES_COMMON) $(SOURCES_$(PARTYPE)) 


OBJECTS = $(SOURCES:.cpp=.o)

INCLUDES = -I$(CURDIR) -I$(CURDIR)/External $(CPPUNIT_INCLUDES)


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
CXXFLAGS = $(CXXFLAGS_$(PARTYPE)_$(TARGET)) $(INCLUDES) -DMPICH_IGNORE_CXX_SEEK 
LINKFLAGS = $(LINKFLAGS_$(PARTYPE)_$(TARGET)) $(CPPUNIT_LINK_ARGS) -lrt 

ifeq ($(DEPS),1)
  DEPCOMP = depmode=$(DEPMODE) source=$< object=$@ $(MAKEFILEDIR)depcomp
else
  DEPCOMP =
endif

ifeq ($(LCOV),1)
  ifeq ($(TARGET),DEBUG)
    $(info # Creating LCOV) 
    CXXFLAGS += -fno-default-inline -fno-inline -fprofile-arcs -ftest-coverage
    LINKFLAGS += -fprofile-arcs -ftest-coverage
  endif
endif

VTK_MAKEFILE := $(realpath $(MAKEFILEDIR)vtkwriter.mk ) 
include $(VTK_MAKEFILE)

$(BINARY): .deps $(OBJECTS) cppunit
	@echo "building $(BINARY)"
	$(CXX) $(LINKFLAGS) -o $@ $(OBJECTS)
	$(shell if [ ! -f $(BINARYbasename) -o -L $(BINARYbasename) ]; then ln -sf $(BINARY) $(BINARYbasename); fi)

%.o: %.cpp
	$(DEPCOMP) $(CXX) $(CXXFLAGS) -c $< -o $@

.PHONY: help all clean cleanall depend cfg_list cfg_help test lcov

help:
	@echo "Current configuration: $(CFG_INFO)"
	@echo
	@echo "options:"
	@echo "make CFG=...                 choose compiler (see make cfg_list)"
	@echo "make PARTYPE=SEQ | PAR       sequential or parallel (MPI) build"
	@echo "make TARGET=DEBUG | RELEASE  debugging / non optimized build"
	@echo "make DEPS=0 | 1              compile using dependancy tracking"
	@echo "make TESTS=0 | 1             compile with or without unit tests (and build cppunit if neccessary)"
	@echo "make VTK=0 | 1	            compile with or without the vtk output plugin" 
	@echo "make LCOV=0 | 1              compile with support for lcov (code-coverage)"
	@echo
	@echo "targets:"
	@echo "make all         build $(BINARY)"
	@echo "make cfg_list    list available configurations"
	@echo "make cfg_help    show help for the selected configuration"
	@echo "make clean       delete object files"
	@echo "make cleanall    delete object files and executable"
	@echo "make depend      generate makefile dependency list"
	@echo "make lcov        create code coverage information"
	@echo "make lcov-clean  reset profile data"

all: $(BINARY)

test: 
	./$(BINARY) -t

#list available configurations
cfg_list:
	@ls -1 ../makefile/cfg | sed -e 's/.in//'
#print info about the selected configuration
cfg_help:
	@echo "Current configuration: $(CFG_INFO)"
	@echo 
	@echo "Available PARTYPEs:"
	@echo -e $(CFG_PARTYPES)
	@echo "Available TARGETs:"
	@echo -e $(CFG_TARGETS)
	

clean:
	find -type f -name '*.o' -delete
	find -path '*/.deps/*' -type f -delete
	find -type d -name .deps -delete
	find -type f -name '*.gcov' -delete
	find -type f -name '*.gcda' -delete
	find -type f -name '*.gcno' -delete

clear: cleanall

cleanall: clean cppunit_clean
	rm -f $(BINARY)
	find -type f -name '*~' -delete
	if [ -L $(BINARYbasename) ]; then rm -f $(BINARYbasename); fi

lcov:
	lcov --directory $(CURDIR) -b $(CURDIR) --capture --output-file Mardyn.info
	genhtml --output-director $(CURDIR)/../lcov/ Mardyn.info

lcov-clean:
	 lcov --directory $(CURDIR) -b $(CURDIR) --zerocounters

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

