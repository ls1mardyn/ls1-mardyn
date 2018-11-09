# lz4 library
option(ENABLE_LZ4 "Use lz4 library" OFF)
if(ENABLE_LZ4)
    message(STATUS "Using LZ4.")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DMARDYN_LZ4")

    # Enable ExternalProject CMake module
    include(ExternalProject)

    # Download and install lz4
    ExternalProject_Add(
            lz4
            GIT_REPOSITORY https://github.com/lz4/lz4.git
            GIT_TAG origin/dev
            #URL https://github.com/LZ4/LZ4/archive/0c3d8b07a2e38940057fafd21b98645cb074e729.zip # zip option
            #${CMAKE_SOURCE_DIR}/libs/googletest-master.zip # bundled option
            #URL_HASH MD5=6e70656897167140c1221eecc6ad872d
            BUILD_BYPRODUCTS lz4/liblz4.so
            PREFIX contrib/cmake_unofficial
            # Disable install step
            INSTALL_COMMAND ""
            CMAKE_ARGS
            -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
            -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
        #     -DBUILD_TESTS=OFF
        #     -DBUILD_EXAMPLES=OFF
        #     -DENABLE_ADDRESS_SANITIZER=${ENABLE_ADDRESS_SANITIZER}
    )

    # Get lz4 source and binary directories from CMake project
    ExternalProject_Get_Property(lz4 source_dir binary_dir)

    # Create a liblz4 target to be used as a dependency by the program
    add_library(liblz4 IMPORTED SHARED GLOBAL)
    add_dependencies(liblz4 lz4)

    set_target_properties(liblz4 PROPERTIES
            "IMPORTED_LOCATION" "${binary_dir}/src/lz4/liblz4.a"
            "IMPORTED_LINK_INTERFACE_LIBRARIES" "${CMAKE_THREAD_LIBS_INIT}"
            )

    # I couldn't make it work with INTERFACE_INCLUDE_DIRECTORIES
    include_directories(SYSTEM
            "${source_dir}/src"
            "${source_dir}/libs/spdlog-0.16.3/include"
            )

    set(lz4_lib "liblz4")
else()
    message(STATUS "Not using LZ4.")
    set(lz4_lib "")
endif()
