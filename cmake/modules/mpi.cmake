# mpi
option(ENABLE_MPI "Enable MPI" ON)
if(ENABLE_MPI)
    message(STATUS "MPI Enabled")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DENABLE_MPI")
else()
    message(STATUS "MPI Disabled")
endif()