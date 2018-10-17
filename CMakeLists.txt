# disable in source builds
set(CMAKE_DISABLE_SOURCE_CHANGES ON)
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)

project(MarDyn)

cmake_minimum_required(VERSION 3.6)

set(CMAKE_CXX_STANDARD 14)

# cmake < 3.6 ignores "set(CMAKE_CXX_STANDARD 11)" when using icpc
if (CMAKE_VERSION VERSION_LESS "3.6" AND CMAKE_CXX_COMPILER_ID MATCHES "Intel")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
endif ()

set(CMAKE_INCLUDE_CURRENT_DIR ON)

# let ccmake and cmake-gui offer the default build type options
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug;Release;RelWithDebInfo;MinSizeRel")

# set Release as the default build type if it is not yet set.
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING
            "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif(NOT CMAKE_BUILD_TYPE)

# the -dynamic flag seems to be needed on hazelhen. no idea what it does, but it fixes some GLIBC linker issue
# this was done following stackoverflow.com/questions/28198965/how-to-resolve-tls-get-addrglibc-2-3-error-on-compiling-gcc-itself

if(CMAKE_CXX_COMPILER_ID MATCHES Cray)
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -dynamic")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -dynamic")
    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -dynamic")
    set(CMAKE_CXX_FLAGS_MINSIZEREL "${CMAKE_CXX_FLAGS_MINSIZEREL} -dynamic")
endif()

# generates test target
ENABLE_TESTING()

# modules
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/modules/)

# ----- first: include option files
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

# ----- second: include libraries, etc.
# include armadillo module
include(armadillo)
# include vtk module
include(vtk)
# include cppunit module
include(cppunit)
# include bhfmm module
include(bhfmm)

# ----- third: include autopas
include(autopas)


# add mardyn
ADD_SUBDIRECTORY(src)

# add libraries
ADD_SUBDIRECTORY(libs)
