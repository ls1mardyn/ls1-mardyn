

# ---- PRECISION ----
# list of available options
set(PRECISION_OPTIONS "DOUBLE;SINGLE;MIXED")
# set instruction set type
set(PRECISION "DOUBLE" CACHE STRING "Precision to use (${PRECISION_OPTIONS}).")
# let ccmake and cmake-gui offer the options
set_property(CACHE PRECISION PROPERTY STRINGS ${PRECISION_OPTIONS})
if (PRECISION MATCHES "^DOUBLE$")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DMARDYN_DPDP")
elseif(PRECISION MATCHES "^SINGLE$")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DMARDYN_SPSP")
elseif(PRECISION MATCHES "^MIXED$")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DMARDYN_SPDP")
else()
    message(FATAL_ERROR "wrong precision option ")
endif()

# ---- PROFILING ----
option(ENABLE_GPROF "Use the GNU profiler gprof (Only supported by GNU/Clang compiler)" OFF)
if(ENABLE_GPROF)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
    if(NOT (CMAKE_CXX_COMPILER_ID MATCHES "Clang" OR CMAKE_CXX_COMPILER_ID MATCHES "GNU"))
        message(WARNING "${CMAKE_CXX_COMPILER_ID} may not support the -pg option.\n(Only use GNU/Clang with the GNU profiler!)\n")
    endif()
endif()

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DQUICKSCHED")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DTASKTIMINGPROFILE")


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DTIMERS")

if (CMAKE_CXX_COMPILER_ID MATCHES "Intel")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fp-model precise")
endif()

option(DEBUG_DECOMP "Activates verbose debugging for the decompositions" OFF)
if(DEBUG_DECOMP)
    add_definitions(-DDEBUG_DECOMP)
endif()

option(REDUCED_MEMORY_MODE "Activates the reduced memory mode" OFF)
if(REDUCED_MEMORY_MODE)
    add_definitions(-DENABLE_REDUCED_MEMORY_MODE=1)
endif()