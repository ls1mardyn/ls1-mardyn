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
                URL ${MARDYN_SOURCE_DIR}/libs/ALL_d16796dc.zip
                URL_HASH MD5=956315034d9d46d4d03741a6669a6dec
        )
    else()
        set(ALLRepoPath https://gitlab.version.fz-juelich.de/SLMS/loadbalancing.git)
        if (GIT_SUBMODULES_SSH)
            set(ALLRepoPath git@gitlab.version.fz-juelich.de:10022/SLMS/loadbalancing.git)
        endif ()

        FetchContent_Declare(
                allfetch
                GIT_REPOSITORY ${ALLRepoPath}
                GIT_TAG d16796dc
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