

# ---- PRECISION ----
# list of available options
set(PRECISION_OPTIONS "DOUBLE;SINGLE;MIXED")
# set instruction set type
set(PRECISION "DOUBLE" CACHE STRING "Precision to use (${PRECISION_OPTIONS}).")
# let ccmake and cmake-gui offer the options
set_property(CACHE PRECISION PROPERTY STRINGS ${PRECISION_OPTIONS})
if (PRECISION MATCHES "^DOUBLE$")
    message(STATUS "Using double precision.")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DMARDYN_DPDP")
elseif(PRECISION MATCHES "^SINGLE$")
    message(STATUS "Using single precision.")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DMARDYN_SPSP")
elseif(PRECISION MATCHES "^MIXED$")
    message(STATUS "Using mixed precision.")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DMARDYN_SPDP")
else()
    message(FATAL_ERROR "wrong precision option ")
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
    message(STATUS "Using reduced memory mode.")
    add_definitions(-DENABLE_REDUCED_MEMORY_MODE=1)
endif()