
$(info included quicksched.mk!)

#TODO remove when both quicksched versions are supported or introduce check for pthreads
ifneq ($(OPENMP),1)
  $(error ERROR: Quicksched needs OpenMP. Please run again with OPENMP=1)
endif

QUICKSCHED_SOURCES= $(shell find ../dependencies-external/quicksched -name '*.c')

ifeq ($(QUICKSCHED_SOURCES),)
  $(error ERROR: Quicksched source files not found in dependencies-external/quicksched !)
endif

#TODO check for custom paths in variables
INCLUDES += -I../dependencies-external/quicksched/src
OBJECTS_C = $(QUICKSCHED_SOURCES:.c=.o)
OBJECTS += $(OBJECTS_C)

#TODO add switches for DTIMERS maybe also PRINT_SCHEDULING_TIMINGS
CXXFLAGS += -DQUICKSCHED
#CXXFLAGS += -DTIMERS

$(info setting C standard to gnu11!)
CFLAGS += -std=gnu11
