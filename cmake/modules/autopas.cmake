# autopas library
option(ENABLE_AUTOPAS "Use autopas library" OFF)
if (ENABLE_AUTOPAS)
    message(STATUS "Using AutoPas.")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DMARDYN_AUTOPAS")

    # Enable ExternalProject CMake module
    include(FetchContent)

    # Select https (default) or ssh path.
    set(autopasRepoPath https://github.com/AutoPas/AutoPas.git)
    if (GIT_SUBMODULES_SSH)
        set(autopasRepoPath git@github.com:AutoPas/AutoPas.git)
    endif ()

    set(AUTOPAS_BUILD_TESTS OFF CACHE BOOL "" FORCE)
    set(AUTOPAS_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
    set(AUTOPAS_ENABLE_ADDRESS_SANITIZER ${ENABLE_ADDRESS_SANITIZER} CACHE BOOL "" FORCE)
    set(AUTOPAS_OPENMP ${OPENMP} CACHE BOOL "" FORCE)
    set(spdlog_ForceBundled ON CACHE BOOL "" FORCE)
    set(Eigen3_ForceBundled ON CACHE BOOL "" FORCE)
    # translate ls1 vectorization settings to AutoPas'
    set(AUTOPAS_USE_VECTORIZATION ${USE_VECTORIZATION} CACHE BOOL "Set via USE_VECTORIZATION" FORCE)
    set(AUTOPAS_VECTOR_INSTRUCTIONS ${VECTOR_INSTRUCTIONS} CACHE STRING "Set via VECTOR_INSTRUCTIONS_OPTIONS" FORCE)

    # Merge of fixSpacialLocksResize 2023-05-08; This includes the rewrites of remainder traversal and iterators.
    set(AUTOPAS_TAG e47c9257bb8d0c6ab680d02823b1e8749b25369a CACHE STRING "AutoPas Git tag or commit id to use")

    FetchContent_Declare(
            autopasfetch
            GIT_REPOSITORY ${autopasRepoPath}
            GIT_TAG ${AUTOPAS_TAG}
    )

    # Get autopas source and binary directories from CMake project
    FetchContent_GetProperties(autopasfetch)

    if (NOT autopasfetch_POPULATED)
        FetchContent_Populate(autopasfetch)

        add_subdirectory(${autopasfetch_SOURCE_DIR} ${autopasfetch_BINARY_DIR} EXCLUDE_FROM_ALL)
    endif ()

    set(AUTOPAS_LIB "autopas")
else ()
    message(STATUS "Not using AutoPas.")
    set(AUTOPAS_LIB "")
endif ()
