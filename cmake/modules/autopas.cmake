# autopas library
option(ENABLE_AUTOPAS "Use autopas library" OFF)
if(ENABLE_AUTOPAS)
    message(STATUS "Using AutoPas.")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DMARDYN_AUTOPAS")

    # Enable ExternalProject CMake module
    include(ExternalProject)

    # Download and install autopas
    ExternalProject_Add(
            autopas
            GIT_REPOSITORY https://github.com/AutoPas/AutoPas.git
            #URL https://github.com/AutoPas/AutoPas/archive/0c3d8b07a2e38940057fafd21b98645cb074e729.zip # zip option
            #${CMAKE_SOURCE_DIR}/libs/googletest-master.zip # bundled option
            #URL_HASH MD5=6e70656897167140c1221eecc6ad872d
            BUILD_BYPRODUCTS ${CMAKE_CURRENT_BINARY_DIR}/autopas/src/autopas/libautopas.a
            PREFIX ${CMAKE_CURRENT_BINARY_DIR}/autopas
            # Disable install step
            INSTALL_COMMAND ""
            CMAKE_ARGS
            -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
            -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
            -DBUILD_TESTS=OFF
            -DBUILD_EXAMPLES=OFF
    )

    # Get autopas source and binary directories from CMake project
    ExternalProject_Get_Property(autopas source_dir binary_dir)

    # Create a libautopas target to be used as a dependency by the program
    add_library(libautopas IMPORTED STATIC GLOBAL)
    add_dependencies(libautopas autopas)

    set_target_properties(libautopas PROPERTIES
            "IMPORTED_LOCATION" "${binary_dir}/src/autopas/libautopas.a"
            "IMPORTED_LINK_INTERFACE_LIBRARIES" "${CMAKE_THREAD_LIBS_INIT}"
            )

    # I couldn't make it work with INTERFACE_INCLUDE_DIRECTORIES
    include_directories(SYSTEM
            "${source_dir}/src"
            "${source_dir}/libs/spdlog-0.16.3/include"
            )

    set(autopas_lib "libautopas")
else()
    message(STATUS "Not using AutoPas.")
    set(autopas_lib "")
endif()