# lz4 library
if(ENABLE_LZ4)
    message(STATUS "Using LZ4.")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DENABLE_LZ4")

    # Enable ExternalProject CMake module
    include(ExternalProject)

    set(LZ4_SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/lz4)
    set(LZ4_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/libs/lz4)
    # Download and install lz4
    ExternalProject_Add(
            lz4
            GIT_REPOSITORY https://github.com/lz4/lz4.git
            GIT_TAG dev
            SOURCE_DIR ${LZ4_SOURCE_DIR}
            BINARY_DIR ${LZ4_BINARY_DIR}
            CONFIGURE_COMMAND cmake ${LZ4_SOURCE_DIR}/contrib/cmake_unofficial -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER} -DBUILD_STATIC_LIBS=ON
            BUILD_COMMAND make
            INSTALL_COMMAND ""    
    )
    # Create a liblz4 target to be used as a dependency by the program
    add_library(liblz4 IMPORTED STATIC)
    add_dependencies(liblz4 lz4)

    # set include directory associated with this target
    include_directories(${LZ4_SOURCE_DIR}/lib)
    set(LZ4_LIB ${LZ4_BINARY_DIR}/liblz4.a)
else()
    message(STATUS "Not using LZ4.")
    set(LZ4_LIB "")
endif()
