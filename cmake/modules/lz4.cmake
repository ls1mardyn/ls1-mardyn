# lz4 library
if (NOT ${LZ4_ENABLED})
    set(LZ4_ENABLED OFF)
endif()
option(ENABLE_LZ4 "Use lz4 library" ${LZ4_ENABLED})
if(ENABLE_LZ4)
    message(STATUS "Using LZ4.")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DMARDYN_LZ4")

    # Enable ExternalProject CMake module
    include(ExternalProject)

    set(LZ4_SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/lz4)
    set(LZ4_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/lib)
    # Download and install lz4
    ExternalProject_Add(
            lz4
            GIT_REPOSITORY https://github.com/lz4/lz4.git
            GIT_TAG dev
            SOURCE_DIR ${LZ4_SOURCE_DIR}
            BINARY_DIR ${LZ4_BINARY_DIR}
            INSTALL_COMMAND ""
            CMAKE_ARGS ${LZ4_SOURCE_DIR}/contrib/cmake_unofficial
            -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
    )

    # Get lz4 source and binary directories from CMake project
    ExternalProject_Get_Property(lz4 source_dir binary_dir)

    # Create a liblz4 target to be used as a dependency by the program
    add_library(liblz4 IMPORTED SHARED GLOBAL)
    add_dependencies(liblz4 lz4)

#     set_target_properties(liblz4 PROPERTIES
#         INTERFACE_INCLUDE_DIRECTORIES ${LZ4_SOURCE_DIR}/lib
#     )
    include_directories(
        ${LZ4_SOURCE_DIR}/lib
    )
    set(LZ4_LIB ${LZ4_BINARY_DIR}/liblz4.so)
else()
    message(STATUS "Not using LZ4.")
    set(LZ4_LIB "")
endif()
