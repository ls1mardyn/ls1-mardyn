# persistent communication
option(ENABLE_PERSISTENT "Enable Persistent" OFF)
if(ENABLE_PERSISTENT)
    message(STATUS "Persistent Communication Enabled")
    # possibly add a check for a MPI version of at least 4.0
    find_package(MPI REQUIRED)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DENABLE_PERSISTENT")
else()
    message(STATUS "Persistent Disabled")
endif()