#
# On windows, store location of Visual Studio compiler
# into the environment. This will be picked up by nvcc,
# even without explicitly being passed.
# On Linux, use whatever gcc is in the current path
# (so leave compiler-bindir undefined):
#
NVCC = nvcc
RM = rm

#
# Similar for OPENCC_FLAGS and PTXAS_FLAGS.
# These are simply passed via the environment:
#
export OPENCC_FLAGS :=
export PTXAS_FLAGS := -fastimul

NVCCFLAGS = -arch=sm_20

#
# Define the application;
# for each object file, there must be a
# corresponding .c or .cpp or .cu file:
#
SOURCES = $(wildcard *.cu) $(wildcard *.cpp) $(wildcard *.c)
OBJECTS = $(addsuffix .o,$(basename $(SOURCES)))
APP = cuda_tests

$(APP) : $(OBJECTS)
	$(NVCC) $(OBJECTS) $(LDFLAGS) -o $@

#
# cuda and C/C++ compilation rules, with
# dependency generation:de	
#
%.o : %.cpp
	$(NVCC) -c $< $(NVCCFLAGS) -o $@
	$(NVCC) -M $< $(NVCCFLAGS) > $@.dep
%.o : %.c
	$(NVCC) -c $< $(NVCCFLAGS) -o $@
	$(NVCC) -M $< $(NVCCFLAGS) > $@.dep
%.o : %.cu
	$(NVCC) -c $< $(NVCCFLAGS) -o $@
	$(NVCC) -ptx $< $(NVCCFLAGS)
	$(NVCC) -M $< $(NVCCFLAGS) > $@.dep
#
# Pick up generated dependency files, and
# add /dev/null because gmake does not consider
# an empty list to be a list:
#
include $(wildcard *.dep) /dev/null

#
# Cleanup:
#
clean :
	$(RM) $(OBJECTS) *.dep