# Find out the base directory
MAKEFILEDIR := $(dir $(realpath $(lastword $(MAKEFILE_LIST))))

.PHONY: all clean clean_deps cleanall default help test
default: all

# Check current directory
CURDIR ?= $(realpath $(PWD) )
SRCDIR := $(CURDIR)
ifeq ($(findstring src,$(SRCDIR)),)
  $(error ERROR: make has to be called from the ls1-MarDyn source directory, current directory is $(CURDIR))
endif

# compiler configuration
# Show list of all available configurations with make cfg_list
CFG ?= gcc

# parallelization type
# Show list of all available parallelization types for the current config with cfg_help
# Each config should at least provide the two parallelization types SEQ and PAR
PARTYPE ?= PAR

# choose target
# Show list of all available targets for the current config with cfg_help
# Each config should at least provide the two targets RELEASE and DEBUG.
TARGET ?= RELEASE

# compile using depmode
# create dependency information for faster recompilation in case of source file changes
# 0: disbaled, 1: enabled
DEPS ?= 0

# build unit tests
# 0: disbaled, 1: enabled
UNIT_TESTS ?= 0

# build with VTK support
# 0: disbaled, 1: enabled
VTK ?= 0

# build with INSITU support
# 0: disbaled, 1: enabled
INSITU ?= 0

# bind to TBB proxy malloc
# 0: disbaled, 1: enabled
TBB ?= 0

# clang's address sanitizer
ADDRESS_SANITIZER ?= 0

# clang's thread sanitizer
THREAD_SANITIZER ?= 0

# select vectorization mode:
# see "make help_vect"
VECTORIZE_CODE ?= SSE
vec_options := NOVEC SSE AVX AVX2 KNL_MASK KNL_G_S SKX_MASK SKX_G_S SSEAMD AVXAMD
deprecated_vec_options := KNC_MASK KNC_G_S
novec_options := AOS SOA
ifeq ($(VECTORIZE_CODE), $(filter $(VECTORIZE_CODE), $(novec_options)))
  $(warning WARNING: AOS and SOA options are now combined to NOVEC. Specify --legacy-cell-processor to run MarDyn to replicate the behavior of AOS.")
  VECTORIZE_CODE = NOVEC
else ifeq ($(VECTORIZE_CODE), $(filter $(VECTORIZE_CODE), $(deprecated_vec_options)))
  $(error ERROR: deprecated vectorize code option chosen. see "make help_vect")
else ifeq (, $(filter $(VECTORIZE_CODE), $(vec_options)))
  $(error ERROR: invalid vectorize code option chosen. see "make help_vect")
endif

PRECISION ?= DOUBLE
# SINGLE - SPSP - computation carried out in single precision, accumulation carried out in single precision
# MIXED  - SPDP - computation carried out in single precision, accumulation carried out in double precision
# DOUBLE - DPDP - computation carried out in double precision, accumulation carried out in double precision
precision_options := SINGLE MIXED DOUBLE
ifeq (, $(filter $(PRECISION), $(precision_options)))
  $(error ERROR: invalid PRECISION option chosen. Set to SINGLE, MIXED or DOUBLE.)
endif
ifeq ($(PRECISION), MIXED)
  $(info # WARNING)
  $(info # WARNING)
  $(info # WARNING: MIXED precision in development! )
  $(info # WARNING)
  $(info # WARNING)
endif

ifeq ($(PRECISION), SINGLE)
  $(info # WARNING)
  $(info # WARNING)
  $(info # WARNING: you are selecting SINGLE precision instead of DOUBLE. Proceed at own risk! )
  $(info # WARNING)
  $(info # WARNING)
endif

# build ls1-MarDyn for the reduced memory mode
# If enabled the reduced memory mode reduces the amount of memory required by
# ls1-MarDyn for very large simulations.
# However, not all features will work in this mode.
# 0: disbaled, 1: enabled
REDUCED_MEMORY_MODE ?= 0

# allow FFT acceleration of FMM (requires FFTW!) ?
ENABLE_FMM_FFT ?= 0

# allow acceleration of MPI communication through overlapping communication ?
ENABLE_OVERLAPPING ?= 0

# allow acceleration of MPI communication through overlapping communication - experimental
ENABLE_ADVANCED_OVERLAPPING ?= 0

# enable internal backtracing for SIGSEGV
SIGHANDLER ?= 0

# allow OpenMP? Leave off for now, unless you are developing it.
OPENMP ?= 0

PAPI ?= 0

ENABLE_QUICKSCHED ?= 0

# use armadillo?
ARMADILLO ?= 0

ARMADILLO_NOWRAP ?= 1

ARMADILLO_NOLINK ?= 0

ifeq ($(ARMADILLO), 1)
	CXXFLAGS += -DMARDYN_ARMADILLO
	ifeq ($(ARMADILLO_NOWRAP), 1)
		CXXFLAGS += -DARMA_DONT_USE_WRAPPER
		ifeq ($(ARMADILLO_NOLINK), 0)
			LDFLAGS += -lblas -llapack
		endif
	else
		LDFLAGS += -larmadillo
	endif
endif


ifeq ($(ENABLE_FMM_FFT), 1)
OPENMP_SIMD = 1
USE_FFTW = 1
else
OPENMP_SIMD = 0
USE_FFTW = 0
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

SOURCES_COMMON = $(shell find ./ -name "*.cpp" | grep -E -v "(parallel/|/tests/|/vtk/|/fft/|AutoPas)")
SOURCES_SEQ = $(shell find parallel/ -name "*.cpp" | grep "DomainDecompBase\|LoadCalc\|Zonal\|ForceHelper" | grep -E -v "/tests/")
SOURCES_PAR = $(shell find parallel/ -name "*.cpp" | grep -E -v "(/tests/|/vtk/|ALLL)")
SOURCES = $(SOURCES_COMMON) $(SOURCES_$(PARTYPE))


OBJECTS = $(SOURCES:.cpp=.o)
INCLUDES = -I$(SRCDIR) -isystem $(SRCDIR)/../libs/rapidxml -isystem $(SRCDIR)/../libs/armadillo

BINARY_BASENAME = MarDyn
ifneq ("$(wildcard $(SRCDIR)/../.git)","")
  $(info Building from a git version.)
  GIT_VERSION = $(shell git rev-parse --short HEAD | sed -e "s/:/-/")
  ifneq ($(GIT_VERSION),)
    BINARY = $(BINARY_BASENAME)_$(GIT_VERSION).$(PARTYPE)_$(TARGET)_$(VECTORIZE_CODE)
  else
    BINARY = $(BINARY_BASENAME).$(PARTYPE)_$(TARGET)_$(VECTORIZE_CODE)
  endif
else
  BINARY = $(BINARY_BASENAME).$(PARTYPE)_$(TARGET)_$(VECTORIZE_CODE)
endif

CXX = $(CXX_$(PARTYPE))
CXXFLAGS += $(CXXFLAGS_$(PARTYPE)_$(TARGET)) $(INCLUDES) -D MPICH_IGNORE_CXX_SEEK
LDFLAGS += $(LINKFLAGS_$(PARTYPE)_$(TARGET))

# additional setup for unit tests
ifneq ($(UNIT_TESTS),0)
  CPPUNIT_MAKE_INCLUDE=$(realpath $(MAKEFILEDIR)/cppunit.mk )
  include $(CPPUNIT_MAKE_INCLUDE)
endif

# set vectorization intrinsics
CXXFLAGS += $(CXXFLAGS_VECTORIZE)
LDFLAGS += $(LDFLAGS_VECTORIZE)

ifeq ($(PRECISION), SINGLE)
CXXFLAGS += -D MARDYN_SPSP
endif
ifeq ($(PRECISION), MIXED)
CXXFLAGS += -D MARDYN_SPDP
endif
ifeq ($(PRECISION), DOUBLE)
CXXFLAGS += -D MARDYN_DPDP
endif

ifeq ($(OPENMP), 1)
CXXFLAGS += $(FLAGS_OPENMP)
LDFLAGS += $(FLAGS_OPENMP)
else
  ifeq ($(OPENMP_SIMD),1)
  CXXFLAGS += $(FLAGS_OPENMP)
  endif
endif

ifeq ($(ENABLE_FMM_FFT),1)
CXXFLAGS += -D FMM_FFT
endif

ifeq ($(SIGHANDLER),1)
CXXFLAGS += -D ENABLE_SIGHANDLER=1
endif

ifeq ($(ENABLE_OVERLAPPING),1)
CXXFLAGS += -D ENABLE_OVERLAPPING
endif

ifeq ($(ENABLE_ADVANCED_OVERLAPPING),1)
CXXFLAGS += -D ADVANCED_OVERLAPPING
endif

ifeq ($(PAPI),1)
CXXFLAGS += -D WITH_PAPI=1
LDFLAGS += -lpapi
endif

ifeq ($(REDUCED_MEMORY_MODE), 1)
CXXFLAGS += -D ENABLE_REDUCED_MEMORY_MODE=1
endif

ifeq ($(DEPS),1)
# create .deps subdirectories containing dependency information
  $(shell find -type d ! -name .deps -execdir mkdir -p {}/.deps \;)
  DEPCOMP = depmode=$(DEPMODE) source=$< object=$@ $(MAKEFILEDIR)depcomp
else
  DEPCOMP =
endif

ifeq ($(INSITU), 1)
	INSITU_MAKEFILE := $(realpath $(MAKEFILEDIR)/insitu.mk )
  include $(INSITU_MAKEFILE)
endif

ifeq ($(VTK), 1)
  VTK_MAKEFILE := $(realpath $(MAKEFILEDIR)/vtkwriter.mk )
  include $(VTK_MAKEFILE)
endif

ifeq ($(TBB), 1)
  TBB_MAKEFILE := $(realpath $(MAKEFILEDIR)/tbbmalloc.mk )
  include $(TBB_MAKEFILE)
endif

ifeq ($(USE_FFTW), 1)
  FFTW_MAKEFILE := $(realpath $(MAKEFILEDIR)/fft.mk )
  include $(FFTW_MAKEFILE)
endif

ifeq ($(ENABLE_QUICKSCHED), 1)
  QUICKSCHED_MAKEFILE := $(realpath $(MAKEFILEDIR)/quicksched.mk )
  include $(QUICKSCHED_MAKEFILE)
endif

ifeq ($(ENABLE_TASK_TIMING_PROFILER), 1)
  CXXFLAGS += -DTASKTIMINGPROFILE
endif

$(BINARY): $(OBJECTS)
	@echo "building $(BINARY)"
	$(CXX) -o $@ $(OBJECTS) $(LDFLAGS)
	$(shell if [ ! -f $(BINARY_BASENAME) -o -L $(BINARY_BASENAME) ]; then ln -sf $(BINARY) $(BINARY_BASENAME); fi)

# make sure C and C++ files are compiled with the same flags excluding -std
%.o: %.c
	$(DEPCOMP) $(CXX) -x c $(CFLAGS) $(filter-out -std=% -pedantic, $(CXXFLAGS)) -c $< -o $@

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
	@echo "make CFG=$(call underline,gcc), ...				choose compiler (see make cfg_list)"
	@echo "make PARTYPE=$(call underline,SEQ) | PAR				choose parallelization type (see cfg_help)"
	@echo "make TARGET=DEBUG | $(call underline,RELEASE)			choose target (see cfg_help)"
	@echo "make DEPS=$(call underline,0) | 1					compile using dependancy tracking"
	@echo "make UNIT_TESTS=$(call underline,0) | 1				compile with or without unit tests (and build cppunit if neccessary)"
	@echo "make VTK=$(call underline,0) | 1	 				compile with or without the vtk output plugin"
	@echo "make VECTORIZE_CODE=NOVEC | $(call underline,SSE) | AVX | ...	compile with vectorization support (see help_vect)"
	@echo "make OPENMP=$(call underline,0) | 1				build with OpenMP support"
	@echo "make PAPI=$(call underline,0) | 1					compile with PAPI counter support"
	@echo "make ENABLE_FMM_FFT=$(call underline,0) | 1			compile using an FFT accelerated FMM method"
	@echo "make SIGHANDLER=$(call underline,0) | 1				compile with or without backtrace capability in case of SIGSEGV"
	@echo
	@echo "targets:"
	@echo "make all		build $(BINARY)"
	@echo "make cfg_list		list available configurations"
	@echo "make cfg_help		show help for the selected configuration"
	@echo "make clean		delete object files"
	@echo "make clean_deps		delete files generated by the depcomp command"
	@echo "make clean_ext		delete object files of external dependencies"
	@echo "make cleanall		delete object files and executable"
	@echo "make depend		generate makefile dependency list"
	@echo "make lcov		create code coverage information"
	@echo "make lcov_clean		reset profile data"
	@echo "make test		run the unit tests"

help_vect:
	@echo "Vectorization help:"
	@echo ""
	@echo "make VECTORIZE_CODE= NOVEC | $(call underline,SSE) | AVX | AVX2 | KNL_MASK | KNL_G_S | SKX_MASK | SKX_G_S"
	@echo " NOVEC    - use the VectorizedCellProcessor with no intrinsics"
	@echo "            For development purposes only, often slower than the LegacyCellProcessor."
	@echo "            You can always run the LegacyCellProcessor by specifying --legacy-cell-processor."
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
	@echo " KNL_MASK - use the VectorizedCellProcessor with AVX512 (knights landing) intrinsics."
	@echo "            Works only on Xeon Phis of the second generation (knights landing)."
	@echo "            Speed-up not yet tested."
	@echo " KNL_G_S  - use the VectorizedCellProcessor with AVX512 (knights landing) intrinsics."
	@echo "            Works only on Xeon Phis of the second generation (knights landing)."
	@echo "            Speed-up not yet tested."
	@echo " SKX_MASK - use the VectorizedCellProcessor with AVX512F intrinsics."
	@echo "            Uses only AVX512F. Only the flags passed to the compiler assume SKX, but not the used intrinsics code."
	@echo "            Speed-up not yet tested."
	@echo " SKX_G_S  - use the VectorizedCellProcessor with AVX512F intrinsics."
	@echo "            Uses only AVX512F. Only the flags passed to the compiler assume SKX, but not the used intrinsics code."
	@echo "            Speed-up not yet tested."
	@echo "            Please make sure that you are using a reasonably updated compiler!"
	@echo ""
	@echo " To determine, whether your current CPU supports SSE, AVX, AVX2, call, respectively: "
	@echo "  cat /proc/cpuinfo | grep flags | tail -n 1 | grep sse3 "
	@echo "  cat /proc/cpuinfo | grep flags | tail -n 1 | grep avx "
	@echo "  cat /proc/cpuinfo | grep flags | tail -n 1 | grep avx2 "
	@echo " If the return line is non-empty, then the CPU supports the respective vectorization mode."
	@echo " Any further encountered problems are likely due to an outdated or misconfigured compiler."
	@echo ""
	@echo " For the best experience, use AVX2 with either the latest GNU or Intel compiler on Intel processors."


all: $(BINARY)


#list available configurations
cfg_list:
	@ls -1 ../makefile/cfg | sed -e '/\.mk$$/d' -e 's/\.in$$//'
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

clean_ext:
	find ../dependencies-external/* -type f -name '*.o' -delete

cleanall: clean clean_deps clean_ext
	$(RM) $(BINARY)
	find -type f -name '*~' -delete
	if [ -L $(BINARY_BASENAME) ]; then $(RM) $(BINARY_BASENAME); fi

clear: cleanall

lcov:
	lcov --directory $(SRCDIR) -b $(SRCDIR) --capture --output-file Mardyn.info
	genhtml --output-director $(SRCDIR)/../lcov/ Mardyn.info

lcov_clean:
	lcov --directory $(SRCDIR) -b $(SRCDIR) --zerocounters


test: $(BINARY_BASENAME)
	@$(eval FILENAME := $(shell ls -l $(BINARY_BASENAME) | tr -s ' ' | cut -d ' ' -f 11))
	$(info Running tests for executable: $(FILENAME))
	@if test $(findstring PAR, $(FILENAME)); then \
		echo parallel execution ; \
		mpirun -n 2 $(SRCDIR)/$(BINARY_BASENAME) -t -d $(SRCDIR)/../test_input/ ; \
		mpirun -n 1 $(SRCDIR)/$(BINARY_BASENAME) -t -d $(SRCDIR)/../test_input/ ; \
	else \
		echo Serial execution ; \
		$(SRCDIR)/$(BINARY_BASENAME) -t -d $(SRCDIR)/../test_input/ ; \
	fi


# include potentially generated dependency information
-include $(shell find -path '*/.deps/*' -type f)
