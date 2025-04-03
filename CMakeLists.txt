# disable in source builds
set(CMAKE_DISABLE_SOURCE_CHANGES ON)
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)

cmake_minimum_required(VERSION 3.24)

project(MarDyn
    VERSION 1.2.0
    DESCRIPTION "ls1-MarDyn - Large Sca1e MoleculAR DYNamics code"
    HOMEPAGE_URL "https://www.ls1-mardyn.de/"
    LANGUAGES C CXX
    )

set(CMAKE_CXX_STANDARD 17)

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
# determine version info and make it available in c++
include(version)
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

# ----- include mamico
include(mamico)

# ----- add adios2 dependency
include(adios2)

# add mardyn
ADD_SUBDIRECTORY(src)

# add libraries
ADD_SUBDIRECTORY(libs)
