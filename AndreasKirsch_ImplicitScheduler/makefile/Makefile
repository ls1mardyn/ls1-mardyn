# compiler configuration
# show list of all available configurations with make cfg_list
CFG ?= gcc

NVCC ?= /usr/local/cuda/bin/nvcc
PTXAS ?= /usr/local/cuda/bin/ptxas

# machine type: SEQ or PAR
PARTYPE ?= SEQ

# choose target: DEBUG|RELEASE
# some configurations have additional targets, see the respective configuration files for details
TARGET ?= DEBUG

# compile using depmode?
DEPS ?= 1

#build unit tests?
TESTS ?= 0

VTK ?= 0

CUDA_CONFIG ?= CUDA_DOUBLE_WBCP
NUM_WARPS ?= 4

MAX_REGISTER_COUNT ?= 63 

MAX_NUM_COMPONENTS ?= 2
MAX_NUM_LJCENTERS ?= 3
MAX_NUM_DIPOLES ?= 3
MAX_NUM_CHARGES ?= 0
MAX_NUM_QUADRUPOLES ?= 0

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
SOURCES_CUDA = $(shell find ./ -name "*.cu" | grep -v "parallel/" | grep -v "/tests/" | grep -v "/vtk/")
SOURCES_CUDA_MODULES = $(shell find ./ -name "*.cum" | grep -v "parallel/" | grep -v "/tests/" | grep -v "/vtk/")
SOURCES_SEQ = $(shell find parallel/ -name "*.cpp" | grep "Dummy" | grep -v "/tests/" | grep -v "/vtk/")
SOURCES_PAR = $(shell find parallel/ -name "*.cpp" | grep -v "Dummy" | grep -v "/tests/" | grep -v "/vtk/")
SOURCES = $(SOURCES_COMMON) $(SOURCES_$(PARTYPE)) 


OBJECTS = $(SOURCES:.cpp=.o) $(SOURCES_CUDA:.cu=.o)
# ptx files are stored in the output directory
PTXS = $(SOURCES_CUDA:.cu=.ptx)

INCLUDES = -I$(CURDIR) -I$(CURDIR)/External $(CPPUNIT_INCLUDES) -I../dependencies-external/nvidia

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
LINKFLAGS = $(LINKFLAGS_$(PARTYPE)_$(TARGET)) $(CPPUNIT_LINK_ARGS) -lrt -lcuda

ifeq ($(DEPS),1)
  DEPCOMP = depmode=$(DEPMODE) source=$< object=$@ $(MAKEFILEDIR)depcomp 
else
  DEPCOMP =
endif

VTK_MAKEFILE := $(realpath $(MAKEFILEDIR)vtkwriter.mk ) 
include $(VTK_MAKEFILE)

$(BINARY): make_config .deps $(OBJECTS) cppunit
	@echo "building $(BINARY)"
	$(CXX) $(LINKFLAGS) -o $@ $(OBJECTS)
	$(shell if [ ! -f $(BINARYbasename) -o -L $(BINARYbasename) ]; then ln -sf $(BINARY) $(BINARYbasename); fi)

%.o: %.cpp
	$(DEPCOMP) $(CXX) $(CXXFLAGS) -c $< -o $@
	
make_config:
	@echo "#define CONFIG_$(CUDA_CONFIG)" > current_make_config.h
	@echo "#define MAX_NUM_WARPS $(NUM_WARPS)" >> current_make_config.h
	@echo "#define MAX_REGISTER_COUNT $(MAX_REGISTER_COUNT)" >> current_make_config.h
	@echo "#define MAX_NUM_COMPONENTS $(MAX_NUM_COMPONENTS)" >> current_make_config.h
	@echo "#define MAX_NUM_LJCENTERS $(MAX_NUM_LJCENTERS)" >> current_make_config.h
	@echo "#define MAX_NUM_DIPOLES $(MAX_NUM_DIPOLES)" >> current_make_config.h
	@echo "#define MAX_NUM_CHARGES $(MAX_NUM_CHARGES)" >> current_make_config.h
	@echo "#define MAX_NUM_QUADRUPOLES $(MAX_NUM_QUADRUPOLES)" >> current_make_config.h
	
	@echo "created temporary current_make_config.h:"
	@cat current_make_config.h
	
	@if ! diff --brief --new-file current_make_config.h cuda/make_config.h; then \
		echo "updating cuda/make_config.h"; \
		cp current_make_config.h cuda/make_config.h; \
	else \
		echo "no changes in cuda/make_config.h"; \
	fi
	@rm current_make_config.h
	
define CU_BUILDER
$(shell $(NVCC) --generate-dependencies $(2) --output-directory ./cuda | sed -e 's/\.o/\.ptx/g' > $(2:.cu=.deps))
endef

$(call CU_BUILDER,./cuda/kernel.ptx,./cuda/kernel.cu)
include ./cuda/kernel.deps

cuda/kernel.ptx: cuda/kernel.cu

%.ptx: %.cu 
	$(NVCC) -arch=compute_20 --use_fast_math -code=sm_20 --ptx -Xptxas -v $< --output-directory $(<D)
	-$(PTXAS) -arch=sm_20 -maxrregcount $(MAX_REGISTER_COUNT) -v $@
%.o: %.ptx
	ld -r -b binary -o $@ $<

.PHONY: help all clean cleanall depend cfg_list test make_config

help:
	@echo "configuration: $(CFG_INFO)"
	@echo
	@echo "options:"
	@echo "	CFG=...                 choose compiler (see make cfg_list)"
	@echo " PARTYPE=SEQ | PAR       sequential or parallel (MPI) build"
	@echo " TARGET=DEBUG | RELEASE  debugging / non optimized build"
	@echo " DEPS=0 | 1              compile using dependancy tracking"
	@echo " TESTS=0 | 1             compile with or without unit tests (and build cppunit if neccessary)"
	@echo " VTK=0 | 1	            compile with or without the vtk output plugin" 
	@echo
	@echo "targets:"
	@echo "make all         build $(BINARY)"
	@echo "make cfg_list    list available configurations"
	@echo "make clean       delete object files"
	@echo "make cleanall    delete object files and executable"
	@echo "make depend      generate makefile dependency list"

all: $(BINARY)

test: 
	./$(BINARY) -t

#list available configurations
cfg_list:
	@ls -1 ../makefile/cfg | sed -e 's/.in//'

clean:
	find -type f -name '*.o' -delete
	find -path '*/.deps/*' -type f -delete
	find -type d -name .deps -delete
	find -type f -name '*.ptx' -delete

clear: cleanall

cleanall: clean cppunit_clean
	rm -f $(BINARY)
	find -type f -name '*~' -delete
	if [ -L $(BINARYbasename) ]; then rm -f $(BINARYbasename); fi

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
	
