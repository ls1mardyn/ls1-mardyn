# Special compile-rule for the VectorizedCellProcessor to enforce inlining.
# NOTE: depending on the installation settings, mpicxx may not actually use icpc! 
# In that case, an error may occur here!

# C++17
CXXFLAGS += -std=c++17 -fp-model precise


# Vectorization settings:
#########################################
ifeq ($(VECTORIZE_CODE),SSE)
CXXFLAGS_VECTORIZE = -msse3
endif
ifeq ($(VECTORIZE_CODE),AVX)
CXXFLAGS_VECTORIZE = -mavx
endif
ifeq ($(VECTORIZE_CODE),AVX2)
CXXFLAGS_VECTORIZE = -march=core-avx2 -fma
# march=... should enable -fma automatically, but we will ensure it.
endif
ifeq ($(VECTORIZE_CODE),KNL_MASK)
CXXFLAGS_VECTORIZE = -xMIC-AVX512
endif
ifeq ($(VECTORIZE_CODE),KNL_G_S)
CXXFLAGS_VECTORIZE = -xMIC-AVX512 -D__VCP_GATHER__
endif
ifeq ($(VECTORIZE_CODE),SKX_MASK)
CXXFLAGS_VECTORIZE = -xCore-AVX512
endif
ifeq ($(VECTORIZE_CODE),SKX_G_S)
CXXFLAGS_VECTORIZE = -xCore-AVX512 -D__VCP_GATHER__
endif

# OpenMP settings:
#########################################
ifeq ($(OPENMP),1)
FLAGS_OPENMP += -qopenmp
else
  ifeq ($(OPENMP_SIMD),1)
  FLAGS_OPENMP += -qopenmp-simd
  endif
endif

particleContainer/adapter/VectorizedCellProcessor.o: particleContainer/adapter/VectorizedCellProcessor.cpp
	$(DEPCOMP) $(CXX) $(CXXFLAGS) -ipo -inline-forceinline -c $< -o $@
particleContainer/adapter/VCP1CLJWR.o: particleContainer/adapter/VCP1CLJWR.cpp
	$(DEPCOMP) $(CXX) $(CXXFLAGS) -ipo -inline-forceinline -c $< -o $@
#io/CubicGridGeneratorInternal.o: io/CubicGridGeneratorInternal.cpp
#    $(DEPCOMP) $(CXX) $(CXXFLAGS) -ipo -inline-forceinline -c $< -o $@
