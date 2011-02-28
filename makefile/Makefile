# Find out the base directory
MAKEFILEDIR := $(dir $(realpath $(lastword $(MAKEFILE_LIST))))

.PHONY: default
default: all

# Check current directory
CURDIR ?= $(realpath $(PWD) )
ifeq ($(findstring src,$(CURDIR)),)
  $(error ERROR: current directory $(CURDIR) seems not to be a source directory)
endif

# compiler configuration
# Show list of all available configurations with make cfg_list
CFG ?= gcc

# parallelization type
# Show list of all available parallelization types for the current config with cfg_help
# Each config should at least provide the two parallelization types SEQ and PAR
PARTYPE ?= SEQ

# choose target
# Show list of all available targets for the current config with cfg_help
# Each config should at least provide the two targets RELEASE and DEBUG.
TARGET ?= RELEASE

# compile using depmode?
DEPS ?= 0

# build unit tests?
UNIT_TESTS ?= 0

# build with VTK support
VTK ?= 0

#build name and include config file
CFGFILE := $(realpath $(MAKEFILEDIR)cfg/$(CFG).in )

ifneq ($(CFGFILE),)
  $(info # Selected configuration: $(CFG) )
  $(info # Including file $(CFGFILE) )
  include $(CFGFILE)
else
  $(error ERROR: config file for $(CFG) not found - no include )
endif

SOURCES_COMMON = $(shell find ./ -name "*.cpp" | grep -E -v "(parallel/|/tests/|/vtk/)")
SOURCES_SEQ = $(shell find parallel/ -name "*.cpp" | grep "Dummy" | grep -E -v "(/tests/|/vtk/)")
SOURCES_PAR = $(shell find parallel/ -name "*.cpp" | grep -E -v "(/tests/|/vtk/)")
SOURCES = $(SOURCES_COMMON) $(SOURCES_$(PARTYPE)) 


OBJECTS = $(SOURCES:.cpp=.o)
INCLUDES = -I$(CURDIR) -I$(CURDIR)/External

BINARY_BASENAME = MarDyn
ifneq ($(realpath .svn),)
  $(info Building from a svn version.)
  SVN_VERSION = $(shell svnversion | sed -e "s/:/-/")
  ifneq ($(SVN_VERSION),)
    BINARY = $(BINARY_BASENAME)_$(SVN_VERSION).$(PARTYPE)_$(TARGET)
  else
    BINARY = $(BINARY_BASENAME).$(PARTYPE)_$(TARGET)
  endif
else
  BINARY = $(BINARY_BASENAME).$(PARTYPE)_$(TARGET)
endif

CXX = $(CXX_$(PARTYPE))
CXXFLAGS = $(CXXFLAGS_$(PARTYPE)_$(TARGET)) $(INCLUDES) -DMPICH_IGNORE_CXX_SEEK 
LDFLAGS = $(LINKFLAGS_$(PARTYPE)_$(TARGET))

# additional setup for unit tests
ifneq ($(UNIT_TESTS),0)
  CPPUNIT_MAKE_INCLUDE=$(realpath $(MAKEFILEDIR)/cppunit.mk )
  include $(CPPUNIT_MAKE_INCLUDE)
endif


ifeq ($(DEPS),1)
  DEPCOMP = depmode=$(DEPMODE) source=$< object=$@ $(MAKEFILEDIR)depcomp
else
  DEPCOMP =
endif


ifeq ($(VTK), 1)
  VTK_MAKEFILE := $(realpath $(MAKEFILEDIR)/vtkwriter.mk ) 
  include $(VTK_MAKEFILE)
endif

$(BINARY): .deps $(OBJECTS)
	@echo "building $(BINARY)"
	$(CXX) $(LDFLAGS) -o $@ $(OBJECTS)
	$(shell if [ ! -f $(BINARY_BASENAME) -o -L $(BINARY_BASENAME) ]; then ln -sf $(BINARY) $(BINARY_BASENAME); fi)

%.o: %.cpp
	$(DEPCOMP) $(CXX) $(CXXFLAGS) -c $< -o $@

.PHONY: help all clean cleanall depend cfg_list cfg_help lcov

help:
	@echo "Current configuration: $(CFG_INFO)"
	@echo
	@echo "options:"
	@echo "make CFG=...                 choose compiler (see make cfg_list)"
	@echo "make PARTYPE=SEQ | PAR       choose parallelization type (see cfg_help)"
	@echo "make TARGET=DEBUG | RELEASE  choose target (see cfg_help)"
	@echo "make DEPS=0 | 1              compile using dependancy tracking"
	@echo "make UNIT_TESTS=0 | 1        compile with or without unit tests (and build cppunit if neccessary)"
	@echo "make VTK=0 | 1               compile with or without the vtk output plugin" 
	@echo
	@echo "targets:"
	@echo "make all         build $(BINARY)"
	@echo "make cfg_list    list available configurations"
	@echo "make cfg_help    show help for the selected configuration"
	@echo "make clean       delete object files"
	@echo "make cleanall    delete object files and executable"
	@echo "make depend      generate makefile dependency list"
	@echo "make lcov        create code coverage information"
	@echo "make lcov_clean  reset profile data"

all: $(BINARY)


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

cleanall: clean
	rm -f $(BINARY)
	find -type f -name '*~' -delete
	if [ -L $(BINARY_BASENAME) ]; then rm -f $(BINARY_BASENAME); fi

lcov:
	lcov --directory $(CURDIR) -b $(CURDIR) --capture --output-file Mardyn.info
	genhtml --output-director $(CURDIR)/../lcov/ Mardyn.info

lcov_clean:
	 lcov --directory $(CURDIR) -b $(CURDIR) --zerocounters

# create .deps subdirectories containing dependency information
.deps:
ifeq ($(DEPS),1)
	$(shell find -type d ! -name .deps -execdir mkdir -p {}/.deps \;)
endif

# include potentially generated dependency information
-include $(shell find -path '*/.deps/*' -type f)
