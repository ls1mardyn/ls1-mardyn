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
    include_directories(SYSTEM
            "${PROJECT_SOURCE_DIR}/dependencies-external/libxsd"
            "$ENV{XERCES_INCDIR}"
            )
else()
    message(STATUS "VTK Disabled")
endif()
