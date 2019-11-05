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
            GIT_TAG f639d8b77eb62b84ffb3717ca4a3e25f1caaea86
            #GIT_TAG origin/feature/regionParticleIteratorIncrease
            #URL https://github.com/AutoPas/AutoPas/archive/0c3d8b07a2e38940057fafd21b98645cb074e729.zip # zip option
            #${CMAKE_SOURCE_DIR}/libs/googletest-master.zip # bundled option
            #URL_HASH MD5=6e70656897167140c1221eecc6ad872d
            BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/autopas/build
            BUILD_BYPRODUCTS ${CMAKE_CURRENT_BINARY_DIR}/autopas/build/src/autopas/libautopas.a
            PREFIX ${CMAKE_CURRENT_BINARY_DIR}/autopas
            # Disable install step
            INSTALL_COMMAND ""
            CMAKE_ARGS
            -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
            -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
            -DAUTOPAS_BUILD_TESTS=OFF
            -DAUTOPAS_BUILD_EXAMPLES=OFF
            -DAUTOPAS_ENABLE_ADDRESS_SANITIZER=${ENABLE_ADDRESS_SANITIZER}
            -DAUTOPAS_OPENMP=${OPENMP}
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

    # workaround for INTERFACE_INCLUDE_DIRECTORIES requiring existent paths, so we create them here...
    file(MAKE_DIRECTORY ${source_dir}/src)
    file(MAKE_DIRECTORY ${source_dir}/libs/spdlog-1.3.1/include)
    file(MAKE_DIRECTORY ${binary_dir}/libs/eigen-3/include)

    target_include_directories(libautopas SYSTEM INTERFACE
            "${source_dir}/src"
            "${source_dir}/libs/spdlog-1.3.1/include"
            "${binary_dir}/libs/eigen-3/include"
            )

    set(AUTOPAS_LIB "libautopas")
else()
    message(STATUS "Not using AutoPas.")
    set(AUTOPAS_LIB "")
endif()
