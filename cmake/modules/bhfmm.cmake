# fast multipole things
option(ENABLE_FMM_FFT "Enable FFT accelerated FMM" OFF)
if(ENABLE_FMM_FFT)
    message(WARNING "FMM_FFT NOT YET SUPPORTED using cmake")
    #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DVTK")
else()
    message(STATUS "FMM_FFT Disabled")
endif()