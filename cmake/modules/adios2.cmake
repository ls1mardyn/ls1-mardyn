option(ENABLE_ADIOS2 "Enables the ADIOS2 writer type." ON)
option(FIND_PACKAGE_ADIOS2 "Uses find_package to find the adios2 library locally instead of downloading it from github." OFF)

if (ENABLE_ADIOS2)
    message(STATUS "Using Adios2.")

    # No local version -> download from github
    if (NOT FIND_PACKAGE_ADIOS2)
        # Enable ExternalProject CMake module
        include(FetchContent)

        # Select https (default) or ssh path.
        set(adios2RepoPath https://github.com/ornladios/ADIOS2.git)
        if (GIT_SUBMODULES_SSH)
            set(adios2RepoPath git@github.com:ornladios/ADIOS2.git)
        endif ()
        
        set(ADIOS2_BUILD_TESTING OFF CACHE BOOL "" FORCE)
        set(ADIOS2_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
        set(BUILD_TESTING OFF CACHE BOOL "" FORCE)
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
                GIT_TAG "v2.7.1.436"
        )

        # Get adios2 source and binary directories from CMake project
        FetchContent_GetProperties(adios2fetch)

        if (NOT adios2fetch_POPULATED)
            FetchContent_Populate(adios2fetch)
            add_subdirectory(${adios2fetch_SOURCE_DIR} ${adios2fetch_BINARY_DIR})
        endif ()

        set(ADIOS2_LIB "adios2")
        set(ADIOS2_COMPILE_DEFINITION "ENABLE_ADIOS2")

        # custom target to copy bpls to the build directory
        add_custom_target(bplsCopy ALL DEPENDS bpls)
        add_custom_command(TARGET bplsCopy POST_BUILD COMMAND ${CMAKE_COMMAND} -E copy ${adios2fetch_BINARY_DIR}/bin/bpls bpls)

    # use local version
    else ()
        find_package(ADIOS2 REQUIRED)
        if (NOT DEFINED ADIOS2_HAVE_MPI)
            set(ADIOS2_HAVE_MPI OFF)
        endif()
        if ((ADIOS2_HAVE_MPI OR ENABLE_MPI) AND (NOT (ADIOS2_HAVE_MPI AND ENABLE_MPI))) # handish implementation of xor (cmake's EQUAL does not work)
            message(FATAL_ERROR "
            You're using an external ADIOS2.
            ADIOS2_HAVE_MPI set to \"${ADIOS2_HAVE_MPI}\" and ls1 ENALBE_MPI is set to \"${ENABLE_MPI}\".
            For a parallel build, MPI has to be enabled in the packaged ADIOS2 as well as ls1 (ADIOS2_ENABLE_MPI=ON, ENABLE_MPI=ON).
            For a sequential build, deactivate MPI for both.")
        endif()
        set(ADIOS2_LIB "adios2::adios2")
        set(ADIOS2_COMPILE_DEFINITION "ENABLE_ADIOS2")
    endif()
else ()
    message(STATUS "Not using adios2.")
    set(ADIOS2_LIB "")
endif ()
