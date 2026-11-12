# cmake module for adding ALL

option(ENABLE_ALLLBL "Enable ALL load balancing library" OFF)
if(ENABLE_ALLLBL)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DENABLE_ALLLBL")
    message(STATUS "ALL load balancing library support enabled.")

    # Enable ExternalProject CMake module
    include(FetchContent)

    option(ALLLBL_USE_BUNDLED "Use bundled version of ALL load balancing library" ON)
    if(ALLLBL_USE_BUNDLED)
        FetchContent_Declare(
                allfetch
                URL ${MarDyn_SOURCE_DIR}/libs/loadbalancing-v0.9.4.zip
                URL_HASH MD5=0913414f33a0d8b7c6f77ac08e2fbc34
        )
    else()
        set(ALLRepoPath https://gitlab.version.fz-juelich.de/SLMS/loadbalancing.git)
        if (GIT_SUBMODULES_SSH)
            set(ALLRepoPath git@gitlab.version.fz-juelich.de:10022/SLMS/loadbalancing.git)
        endif ()

        FetchContent_Declare(
                allfetch
                GIT_REPOSITORY ${ALLRepoPath}
                GIT_TAG v0.9.4
        )
    endif()

    # Get autopas source and binary directories from CMake project
    FetchContent_GetProperties(allfetch)

    if (NOT allfetch_POPULATED)
        FetchContent_Populate(ALLfetch)

        add_library(ALL INTERFACE)
        target_include_directories(ALL INTERFACE ${allfetch_SOURCE_DIR}/include/)
    endif ()
    set(ALL_LIB "ALL")
else()
    message(STATUS "ALL load balancing library support disabled")
    set(ALL_LIB "")
endif()