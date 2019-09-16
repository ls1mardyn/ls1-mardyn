# vtk things
option(ENABLE_VTK "Enable VTK" OFF)
if(ENABLE_VTK)
    message(STATUS "VTK Enabled")
    find_library(VTK_LIB xerces-c
      HINTS $ENV{XERCES_LIBDIR}
      )

    if(NOT VTK_LIB)
        message(FATAL_ERROR "xerces-c lib not found. Disable VTK support, if you do not need it.")
    endif()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DVTK")

    # if the environment variable pointing to the xerces-c include directory is missing try to deduce it
    if(NOT DEFINED ENV{XERCES_INCDIR})
        string(REGEX REPLACE "/lib/.+" "/include/" VTK_INCDIR ${VTK_LIB})
        if(NOT EXISTS ${VTK_INCDIR})
            message(FATAL_ERROR "Could not find include directory for xerces-c. Please set the environment variable XERCES_INCDIR to it.")
        endif()
    else()
        set(VTK_INCDIR $ENV{XERCES_INCDIR})
    endif()

    include_directories(SYSTEM
            "${PROJECT_SOURCE_DIR}/dependencies-external/libxsd"
            ${VTK_INCDIR}
            )
else()
    message(STATUS "VTK Disabled")
endif()
