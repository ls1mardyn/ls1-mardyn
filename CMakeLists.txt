cmake_minimum_required(VERSION 3.5)
project(MarDyn_trunk)
set(SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)

set(CMAKE_CXX_STANDARD 11)

# custom target basically just uses the makefile
# configure here as needed
add_custom_target(MarDyn COMMAND make -C ${SOURCE_DIR}
        CFG=gcc
        PARTYPE=SEQ
        VECTORIZE_CODE=AVX
        OPENMP=1
        ENABLE_QUICKSCHED=0
        -j12
        UNIT_TESTS=0
        CLION_EXE_DIR=${SOURCE_DIR}
        )

# globs for IDEs to see the necessary code
file(GLOB_RECURSE DEPENDENCIES-EXTERNAL
        "dependencies-external/*.cpp"
        "dependencies-external/*.c"
        "dependencies-external/*.h"
        "dependencies-external/*.hxx"
        )

file(GLOB_RECURSE SOURCE_FILES
        "src/*.cpp"
        "src/*.c"
        "src/*.h"
        "src/*.hpp"
        )

file(GLOB_RECURSE TOOLS
        "tools/*.cpp"
        "src/*.h"
        )


# Manually set flags for an IDE to see what you want it to see
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DENABLE_VECTORIZED_CODE=1")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DMARDYN_DPDP")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -DENABLE_OPENMP")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DQUICKSCHED")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DFMM_FFT")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DTASKTIMINGPROFILE")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUNIT_TESTS")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DENABLE_MPI")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DTIMERS")

#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fopenmp ${MARDYN_FLAGS} ${QUICKSCHED_FLAGS}")
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DTIMERS")

include_directories(${SOURCE_DIR})
include_directories(dependencies-external/quicksched/src)
add_executable(MarDyn_cmake ${DEPENDENCIES-EXTERNAL} ${SOURCE_FILES} ${TOOLS})

