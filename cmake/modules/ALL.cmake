# cmake module for adding ALL

option(ENABLE_ALLLBL "Enable ALL load balancing library" OFF)
if(ENABLE_ALLLBL)
    #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DVTK")
    message(STATUS "ALL load balancing library support enabled")
else()
    message(STATUS "ALL load balancing library support disabled")
endif()