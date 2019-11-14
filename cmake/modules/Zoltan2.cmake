# cmake module for adding Zoltan2

option(ENABLE_ZOLTAN2 "Enable Zoltan2 as load balancing library" OFF)
if(ENABLE_ZOLTAN2)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DENABLE_ZOLTAN2")
    FIND_PACKAGE(Zoltan2 REQUIRED)
    MESSAGE("   Trilinos_INCLUDE_DIRS = ${Trilinos_INCLUDE_DIRS}")
    target_include_directories(trilinos_zoltan2 SYSTEM INTERFACE
            "/usr/include/trilinos"
            )
    message(STATUS "link with \${Zoltan2_LIBRARIES}: ${Zoltan2_LIBRARIES} ")
    message(STATUS "Zoltan2 support enabled")
else()
    message(STATUS "Zoltan2 support disabled")
endif()