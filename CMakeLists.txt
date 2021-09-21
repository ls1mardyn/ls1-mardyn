# disable in source builds
set(CMAKE_DISABLE_SOURCE_CHANGES ON)
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)

cmake_minimum_required(VERSION 3.11)

project(MarDyn)

set(MARDYN_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})
set(MARDYN_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR})

set(CMAKE_CXX_STANDARD 17)
if (NOT CYGWIN)
    set(CMAKE_CXX_EXTENSIONS OFF)
endif ()

set(CMAKE_INCLUDE_CURRENT_DIR ON)

# let ccmake and cmake-gui offer the default build type options
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug;Release;RelWithDebInfo;MinSizeRel")

# set Release as the default build type if it is not yet set.
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING
            "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif(NOT CMAKE_BUILD_TYPE)

# generates test target
ENABLE_TESTING()

# modules
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/modules/)

# ----- first: include option files
option(GIT_SUBMODULES_SSH "Use SSH for git submodules instead of HTTPS" OFF)
# include vectorization module
include(vectorization)
# include OpenMP module
include(OpenMP)
# include options module
include(options)
# include mpi module
include(mpi)
# include module for clang sanitizers
include(clang-sanitizers)
# ccache
include(ccache)

# ----- second: include libraries, etc.
# include armadillo module
include(armadillo)
# include vtk module
include(vtk)
# include cppunit module
include(cppunit)
# include bhfmm module
include(bhfmm)
# include ALL load balanding library
include(ALL)

# ----- third: include autopas
include(autopas)

# ----- resilience plugin dependencies are set here
include(compression)

# add mardyn
ADD_SUBDIRECTORY(src)

# add libraries
ADD_SUBDIRECTORY(libs)

# ARM SVE
include(sve)