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


    option(AUTOPAS_BUILD_TESTS "" OFF)
    option(AUTOPAS_BUILD_EXAMPLES "" OFF)
    option(AUTOPAS_ENABLE_ADDRESS_SANITIZER "" ${ENABLE_ADDRESS_SANITIZER})
    option(AUTOPAS_OPENMP "" ${OPENMP})
    option(spdlog_ForceBundled "" ON)
    option(Eigen3_ForceBundled "" ON)

    FetchContent_Declare(
            autopasfetch
            GIT_REPOSITORY ${autopasRepoPath}
            GIT_TAG 918ce7c356250df60b9a1e217482a0c8d4f7bf0d
    )

    # Get autopas source and binary directories from CMake project
    FetchContent_GetProperties(autopas)

    if (NOT autopasfetch_POPULATED)
        FetchContent_Populate(autopasfetch)

        add_subdirectory(${autopasfetch_SOURCE_DIR} ${autopasfetch_BINARY_DIR} EXCLUDE_FROM_ALL)
    endif ()

    set(AUTOPAS_LIB "autopas")
else ()
    message(STATUS "Not using AutoPas.")
    set(AUTOPAS_LIB "")
endif ()
