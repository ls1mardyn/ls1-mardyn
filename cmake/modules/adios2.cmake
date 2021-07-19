option(ENABLE_ADIOS2 "Enables the ADIOS2 writer type." OFF)
option(FIND_PACKAGE_ADIOS2 "Uses find_package to find the adios2 library." OFF)

if (ENABLE_ADIOS2)
    if (NOT FIND_PACKAGE_ADIOS2)
        message(STATUS "Using Adios2.")

        # Enable ExternalProject CMake module
        include(FetchContent)

        # Select https (default) or ssh path.
        set(adios2RepoPath https://github.com/ornladios/ADIOS2.git)
        if (GIT_SUBMODULES_SSH)
            set(adios2RepoPath git@github.com:ornladios/ADIOS2.git)
        endif ()

        set(ADIOS2_BUILD_TESTING OFF CACHE BOOL "" FORCE)
        set(ADIOS2_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
        set(ADIOS2_USE_Fortran OFF CACHE BOOL "" FORCE)
        set(ADIOS2_USE_Python OFF CACHE BOOL "" FORCE)
        if(NOT ENABLE_MPI)
            set(ADIOS2_USE_MPI OFF CACHE BOOL "" FORCE)
        else()
            set(ADIOS2_USE_MPI ON CACHE BOOL "" FORCE)
        endif()

        FetchContent_Declare(
                adios2fetch
                GIT_REPOSITORY ${adios2RepoPath}
                GIT_TAG "v2.5.0"
        )

        # Get autopas source and binary directories from CMake project
        FetchContent_GetProperties(adios2fetch)

        if (NOT adios2fetch_POPULATED)
            FetchContent_Populate(adios2fetch)
            add_subdirectory(${adios2fetch_SOURCE_DIR} ${adios2fetch_BINARY_DIR})
        endif ()

        set(ADIOS2_LIB "adios2")
    endif()
else ()
    message(STATUS "Not using adios2.")
    set(ADIOS2_LIB "")
endif ()
