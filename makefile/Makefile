# Find out the base directory
MAKEFILEDIR := $(dir $(realpath $(lastword $(MAKEFILE_LIST))))

.PHONY: default help clean cleanall
default: all

# Check current directory
CURDIR ?= $(realpath $(PWD) )
SRCDIR := $(CURDIR)
ifeq ($(findstring src,$(SRCDIR)),)
  $(error ERROR: current directory $(SRCDIR) seems not to be a source directory)
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

SUPPORT_GENERATOR ?= 0

# select vectorization mode:
VECTORIZE_CODE ?= SSE
# AOS      - LegacyCellProcessor
# SSE      - VectorizedCellProcessor with SSE intrinsics (default)
# AVX      - VectorizedCellProcessor with AVX intrinsics
# AVX2     - VectorizedCellProcessor with AVX2 intrinsics
# KNC_MASK - VectorizedCellProcessor with KNC intrinsics and masking
# KNC_G_S  - VectorizedCellProcessor with KNC intrinsics and gather-scatter (not fully optimized yet in this branch)
# SOA      - VectorizedCellProcessor with no intrinsics
# see "make help_vect"
vec_options := AOS SSE AVX AVX2 KNC_MASK KNC_G_S SOA SSEAMD AVXAMD
ifeq (, $(filter $(VECTORIZE_CODE), $(vec_options)))
  $(error ERROR: no valid vectorize code option chosen. see "make help_vect")
endif

# allow FFT acceleration of FMM (requires FFTW!) ?
ENABLE_FMM_FFT ?= 0

# allow acceleration of MPI communication through overlapping communication ?
ENABLE_OVERLAPPING ?= 0

# allow OpenMP? Leave off for now, unless you are developing it.
OPENMP ?= 0

ifeq ($(ENABLE_FMM_FFT), 1)
OPENMP_SIMD = 1
USE_FFTW = 1
endif

#build name and include config file
CFGFILE := $(realpath $(MAKEFILEDIR)cfg/$(CFG).in )

ifneq ($(CFGFILE),)
  $(info # Selected configuration: $(CFG) )
  $(info # Including file $(CFGFILE) )
  include $(CFGFILE)
else
  $(error ERROR: config file for $(CFG) not found - no include )
endif

SOURCES_COMMON = $(shell find ./ -name "*.cpp" | grep -E -v "(parallel/|/tests/|/vtk/|/fft/)")
SOURCES_SEQ = $(shell find parallel/ -name "*.cpp" | grep "Dummy" | grep -E -v "(/tests/|/vtk/)")
SOURCES_SEQ += $(shell find parallel/ -name "*.cpp" | grep "DomainDecompBase" | grep -E -v "/tests/")
SOURCES_PAR = $(shell find parallel/ -name "*.cpp" | grep -E -v "(/tests/|/vtk/)")
SOURCES = $(SOURCES_COMMON) $(SOURCES_$(PARTYPE)) 


OBJECTS = $(SOURCES:.cpp=.o)
INCLUDES = -I$(SRCDIR) -isystem $(SRCDIR)/External

BINARY_BASENAME = MarDyn
ifneq ($(shell realpath .svn),)
  $(info Building from a svn version.)
  SVN_VERSION = $(shell svnversion | sed -e "s/:/-/")
  ifneq ($(SVN_VERSION),)
    BINARY = $(BINARY_BASENAME)_$(SVN_VERSION).$(PARTYPE)_$(TARGET)_$(VECTORIZE_CODE)
  else
    BINARY = $(BINARY_BASENAME).$(PARTYPE)_$(TARGET)_$(VECTORIZE_CODE)
  endif
else
  BINARY = $(BINARY_BASENAME).$(PARTYPE)_$(TARGET)_$(VECTORIZE_CODE)
endif

CXX = $(CXX_$(PARTYPE))
CXXFLAGS += $(CXXFLAGS_$(PARTYPE)_$(TARGET)) $(INCLUDES) -DMPICH_IGNORE_CXX_SEEK -std=c++11
LDFLAGS = $(LINKFLAGS_$(PARTYPE)_$(TARGET)) 

ifeq ($(SUPPORT_GENERATOR), 1)
LDFLAGS += -L../tools/gui/staticlibs -lCubicGridGenerator -lDropletGenerator
CXXFLAGS += -DSUPPORT_GENERATOR
INCLUDES += -I../tools/gui/ -I../tools/gui/src/
endif

# additional setup for unit tests
ifneq ($(UNIT_TESTS),0)
  CPPUNIT_MAKE_INCLUDE=$(realpath $(MAKEFILEDIR)/cppunit.mk )
  include $(CPPUNIT_MAKE_INCLUDE)
endif

ifeq ($(VECTORIZE_CODE),AOS)
CXXFLAGS += -D ENABLE_VECTORIZED_CODE=0
else
CXXFLAGS += -D ENABLE_VECTORIZED_CODE=1
CXXFLAGS += $(CXXFLAGS_VECTORIZE)
LDFLAGS += $(LDFLAGS_VECTORIZE)
endif

ifeq ($(OPENMP), 1)
CXXFLAGS += -DENABLE_OPENMP $(FLAGS_OPENMP)
LDFLAGS += $(FLAGS_OPENMP)
else
  ifeq ($(OPENMP_SIMD),1)
  CXXFLAGS += $(FLAGS_OPENMP)
  endif
endif

ifeq ($(ENABLE_FMM_FFT),1)
CXXFLAGS += -DFMM_FFT
endif

ifeq ($(ENABLE_OVERLAPPING),1)
CXXFLAGS += -DENABLE_OVERLAPPING
endif

ifeq ($(PAPI),1)
CXXFLAGS += -D WITH_PAPI=1
LDFLAGS += -lpapi
endif


ifeq ($(DEPS),1)
  # create .deps subdirectories containing dependency information
  $(shell find -type d ! -name .deps -execdir mkdir -p {}/.deps \;)
  DEPCOMP = depmode=$(DEPMODE) source=$< object=$@ $(MAKEFILEDIR)depcomp
else
  DEPCOMP =
endif


ifeq ($(VTK), 1)
  VTK_MAKEFILE := $(realpath $(MAKEFILEDIR)/vtkwriter.mk ) 
  include $(VTK_MAKEFILE)
endif

ifeq ($(USE_FFTW), 1)
  FFTW_MAKEFILE := $(realpath $(MAKEFILEDIR)/fft.mk ) 
  include $(FFTW_MAKEFILE)
endif

$(BINARY): $(OBJECTS)
	@echo "building $(BINARY)"
	$(CXX) -o $@ $(OBJECTS) $(LDFLAGS)
	$(shell if [ ! -f $(BINARY_BASENAME) -o -L $(BINARY_BASENAME) ]; then ln -sf $(BINARY) $(BINARY_BASENAME); fi)

%.o: %.cpp
	$(DEPCOMP) $(CXX) $(CXXFLAGS) -c $< -o $@

.PHONY: help all clean cleanall depend cfg_list cfg_help lcov clean_deps

hl=$(shell tput smul)
norm=$(shell tput sgr0)
underline=${hl}$(1)${norm}
help:
	@echo "Current configuration: $(CFG_INFO)"
	@echo
	@echo "options (underlined options are default values):"
	@echo "make CFG=$(call underline,gcc), ...                      			choose compiler (see make cfg_list)"
	@echo "make PARTYPE=$(call underline,SEQ) | PAR            			choose parallelization type (see cfg_help)"
	@echo "make TARGET=DEBUG | $(call underline,RELEASE)       			choose target (see cfg_help)"
	@echo "make DEPS=$(call underline,0) | 1                   			compile using dependancy tracking"
	@echo "make UNIT_TESTS=$(call underline,0) | 1             			compile with or without unit tests (and build cppunit if neccessary)"
	@echo "make VTK=$(call underline,0) | 1                    	 		compile with or without the vtk output plugin" 
	@echo "make SUPPORT_GENERATOR=$(call underline,0) | 1   	     		   	compile with support for input generators as dynamic libraries"
	@echo "make VECTORIZE_CODE=AOS | $(call underline,SSE) | AVX | ...		compile with vectorization support (see help_vect)"
	@echo "make PAPI=$(call underline,0) | 1			                   	compile with PAPI counter support"
	@echo "make ENABLE_FMM_FFT=$(call underline,0) | 1				compile using an FFT accelerated FMM method"
	@echo
	@echo "targets:"
	@echo "make all         build $(BINARY)"
	@echo "make cfg_list    list available configurations"
	@echo "make cfg_help    show help for the selected configuration"
	@echo "make clean       delete object files"
	@echo "make clean_deps  delete files generated by the depcomp command"
	@echo "make cleanall    delete object files and executable"
	@echo "make depend      generate makefile dependency list"
	@echo "make lcov        create code coverage information"
	@echo "make lcov_clean  reset profile data"
	@echo "make test  		run the unit tests"

help_vect:
	@echo "Vectorization help:"
	@echo
	@echo "make VECTORIZE_CODE=AOS | SSE | AVX | AVX2 | KNC_MASK | KNC_G_S"
	@echo " AOS      - use the LegacyCellProcessor"
	@echo " SSE      - use the VectorizedCellProcessor with SSE intrinsics"
	@echo "            This is the default option."
	@echo "            Should work on all CPUs newer than 2006 or so."
	@echo "            Speed-up of almost a factor of 2 versus the LegacyCellProcessor."
	@echo " AVX      - use the VectorizedCellProcessor with AVX intrinsics"
	@echo "            Should work on CPUs newer than 2011"
	@echo "            (Intel SandyBridge and AMD Bulldozer are the earliest to support it)."
	@echo "            Speed-up of a factor of 2 to 3 versus the LegacyCellProcessor."
	@echo "            Please make sure that you are using a reasonably updated compiler!"
	@echo "            For AMD Bulldozer, please use the GNU compiler with the additional flags (see cfg/gcc.in)."
	@echo " AVX2     - use the VectorizedCellProcessor with AVX2 intrinsics"
	@echo "            Should work on CPUs newer than 2013"
	@echo "            (Intel Haswell and AMD Excavator are the earliest to support it)."
	@echo "            Speed-up of single digit percentages compared with pure AVX instrinsics."
	@echo "            Please make sure that you are using a reasonably updated compiler!"
	@echo "            For AMD Bulldozer, please use the GNU compiler with the additional flags (see cfg/gcc.in)."
	@echo " KNC_MASK - use the VectorizedCellProcessor with MIC (knights corner) intrinsics."
	@echo "            Works only on Xeon Phis of the first generation (knights corner)."
	@echo "            Speed-up not yet tested."
	@echo " KNC_G_S  - use the VectorizedCellProcessor with MIC (knights corner) intrinsics."
	@echo "            Works only on Xeon Phis of the first generation (knights corner)."
	@echo "            Speed-up not yet tested."
	@echo "            Please make sure that you are using a reasonably updated compiler!"
	@echo " SOA      - use the VectorizedCellProcessor with no intrinsics ("no-vec" mode)"
	@echo "            For development purposes only. Often slower than the LegacyCellProcessor."
	@echo " To determine, whether your current CPU supports SSE, AVX, AVX2, call, respectively: "
	@echo " cat /proc/cpuinfo | grep flags | tail -n 1 | grep sse3 "
	@echo " cat /proc/cpuinfo | grep flags | tail -n 1 | grep avx "
	@echo " cat /proc/cpuinfo | grep flags | tail -n 1 | grep avx2 "
	@echo " If the return line is non-empty, then the CPU supports the respective vectorization mode."
	@echo " Any further encountered problems are likely due to an outdated or misconfigured compiler."
	@echo
	@echo " For the best experience, use AVX2 with either the latest GNU or Intel compiler on Intel processors."


all: $(BINARY)


#list available configurations
cfg_list:
	@ls -1 ../makefile/cfg | sed -e 's/\.in$$//'
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
	find -type f -name '*.gcov' -delete
	find -type f -name '*.gcda' -delete
	find -type f -name '*.gcno' -delete

clean_deps:
	find -path '*/.deps/*' -type f -delete
	find -type d -name .deps -delete

clear: cleanall

cleanall: clean clean_deps
	$(RM) $(BINARY)
	find -type f -name '*~' -delete
	if [ -L $(BINARY_BASENAME) ]; then $(RM) $(BINARY_BASENAME); fi

lcov:
	lcov --directory $(SRCDIR) -b $(SRCDIR) --capture --output-file Mardyn.info
	genhtml --output-director $(SRCDIR)/../lcov/ Mardyn.info

lcov_clean:
	 lcov --directory $(SRCDIR) -b $(SRCDIR) --zerocounters

test:
	$(SRCDIR)/MarDyn -t -d $(SRCDIR)/../test_input/

# include potentially generated dependency information
-include $(shell find -path '*/.deps/*' -type f)
