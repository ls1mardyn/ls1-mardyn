10/11/2011

How to build MarDyn:

make -C src CUDA_CONFIG=[see below} NUM_WARPS=[see below]

CUDA_CONFIG can be eg CUDA_DOUBLE_REFERENCE, CUDA_DOUBLE, CUDA_FLOAT, CUDA_FLOAT_WBCP, CUDA_DOUBLE_WBCP or NO_CUDA
see src/cuda/config.h for all possible configurations

comment out BENCHMARKING in src/cuda/config.h to prevent MarDyn from outputing benchmark data and run normally.

NUM_WARPS=1,8,16
depending on the config

8 or 16 for *_WBCP and high cell densities
1 for other configs and low cell densities
