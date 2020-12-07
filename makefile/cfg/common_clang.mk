CXXFLAGS += -std=c++17
CXXFLAGS_WARNINGS =  -pedantic -Wall -Wextra -Wcast-align -Wcast-qual -Wctor-dtor-privacy -Wdisabled-optimization -Wno-pointer-arith -Wformat=2 -Winit-self -Wmissing-include-dirs -Woverloaded-virtual -Wredundant-decls -Wshadow -Wstrict-overflow -Wundef -Wno-unused # -Wsign-conversion -Wswitch-default -Wold-style-cast -Wsign-promo -Wmissing-declarations 


# address sanitizer:
ifeq ($(ADDRESS_SANITIZER),1)
CXXFLAGS += -fsanitize=address
LDFLAGS += -fsanitize=address
endif

# thread sanitizer:
ifeq ($(THREAD_SANITIZER),1)
CXXFLAGS += -fsanitize=thread
LDFLAGS += -fsanitize=thread
endif


# Vectorization settings:
#########################################
ifeq ($(VECTORIZE_CODE),SSE)
CXXFLAGS_VECTORIZE = -msse3
endif
ifeq ($(VECTORIZE_CODE),AVX)
CXXFLAGS_VECTORIZE = -mavx
endif
ifeq ($(VECTORIZE_CODE),AVX2)
CXXFLAGS_VECTORIZE = -mavx2 -mfma
endif
ifeq ($(VECTORIZE_CODE),KNL_MASK)
CXXFLAGS_VECTORIZE = -march=knl
endif
ifeq ($(VECTORIZE_CODE),KNL_G_S)
CXXFLAGS_VECTORIZE = -march=knl -D__VCP_GATHER__
endif
ifeq ($(VECTORIZE_CODE),KNL_MASK)
CXXFLAGS_VECTORIZE = -march=skylake-avx512
endif
ifeq ($(VECTORIZE_CODE),KNL_G_S)
CXXFLAGS_VECTORIZE = -march=skylake-avx512 -D__VCP_GATHER__
endif

# On AMD BUlldozer:
ifeq ($(VECTORIZE_CODE),SSEAMD)
CXXFLAGS_VECTORIZE = -msse3 -mfma4 -march=bdver1
endif
ifeq ($(VECTORIZE_CODE),AVXAMD)
CXXFLAGS_VECTORIZE = -mavx -mfma4 -march=bdver1
endif

# OpenMP settings:
#########################################
ifeq ($(OPENMP),1)
FLAGS_OPENMP += -fopenmp
else
  ifeq ($(OPENMP_SIMD),1)
  FLAGS_OPENMP += -fopenmp-simd
  endif
endif
