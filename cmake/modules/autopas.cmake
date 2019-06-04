# autopas library
option(ENABLE_AUTOPAS "Use autopas library" OFF)
if(ENABLE_AUTOPAS)
    message(STATUS "Using AutoPas.")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DMARDYN_AUTOPAS")

    # Enable ExternalProject CMake module
    include(ExternalProject)

    # Select https (default) or ssh path.
    set(autopasRepoPath https://github.com/AutoPas/AutoPas.git)
    if(GIT_SUBMODULES_SSH)
        set(autopasRepoPath git@github.com:AutoPas/AutoPas.git)
    endif()

    # Download and install autopas
    ExternalProject_Add(
            autopas
            GIT_REPOSITORY ${autopasRepoPath}
            GIT_TAG origin/master
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
            -DENABLE_ADDRESS_SANITIZER=${ENABLE_ADDRESS_SANITIZER}
            -DOPENMP=${OPENMP}
            -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
    )

    # Get autopas source and binary directories from CMake project
    ExternalProject_Get_Property(autopas source_dir binary_dir)

    # Create a libautopas target to be used as a dependency by the program
    add_library(libautopas IMPORTED STATIC GLOBAL)
    add_dependencies(libautopas autopas)

    # Using target_compile_definitions for imported targets is only possible starting with cmake 3.11, so we use add_definitions here.
    if(OPENMP)
        add_definitions(-DAUTOPAS_OPENMP)
    endif(OPENMP)

    set_target_properties(libautopas PROPERTIES
            "IMPORTED_LOCATION" "${binary_dir}/src/autopas/libautopas.a"
            "IMPORTED_LINK_INTERFACE_LIBRARIES" "${CMAKE_THREAD_LIBS_INIT}"
            )

    # Using INTERFACE_INCLUDE_DIRECTORIES is only possible starting with cmake 3.11, so we use include_directories here!
    include_directories(SYSTEM
            "${source_dir}/src"
            "${source_dir}/libs/spdlog-1.3.1/include"
            )

    set(autopas_lib "libautopas")
else()
    message(STATUS "Not using AutoPas.")
    set(autopas_lib "")
endif()
