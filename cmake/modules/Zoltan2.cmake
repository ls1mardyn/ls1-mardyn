# cmake module for adding Zoltan2

function(check_defined configFile name_of_define RESULT_NAME)
    execute_process(
         COMMAND grep -q -i "#define ${name_of_define}" ${configFile}
         RESULT_VARIABLE define_found
         OUTPUT_QUIET
         )
    set (${RESULT_NAME} ${define_found} PARENT_SCOPE)
endfunction()

option(ENABLE_ZOLTAN2 "Enable Zoltan2 as load balancing library" OFF)
if(ENABLE_ZOLTAN2)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DENABLE_ZOLTAN2")
    FIND_PACKAGE(Zoltan2 REQUIRED)
    MESSAGE(STATUS "Zoltan2_INCLUDE_DIRS = ${Zoltan2_INCLUDE_DIRS}")

    list (GET Zoltan2_LIBRARIES 0 Zoltan2LibraryTarget)

    MESSAGE(STATUS "using Zoltan2_INCLUDE_DIRS")
    target_include_directories(${Zoltan2LibraryTarget} SYSTEM INTERFACE
            "${Zoltan2_INCLUDE_DIRS}"
            )

    # check if we need to include scotch
    check_defined(${Zoltan2_INCLUDE_DIRS}/Zoltan2_config.h HAVE_ZOLTAN2_SCOTCH NEEDS_SCOTCH)
    if(NOT ${NEEDS_SCOTCH})
        message(STATUS "we need scotch")
        target_include_directories(${Zoltan2LibraryTarget} SYSTEM INTERFACE
            "/usr/include/scotch"
        )
    else()
        message(STATUS "skipping scotch, as zoltan2 doesn't need it.")
    endif()

    # check if we need to include parmetis
    check_defined(${Zoltan2_INCLUDE_DIRS}/Zoltan2_config.h HAVE_ZOLTAN2_PARMETIS NEEDS_PARMETIS)
    if(NOT ${NEEDS_PARMETIS})
        message(STATUS "we need parmetis")
        target_include_directories(${Zoltan2LibraryTarget} SYSTEM INTERFACE
                "$ENV{PARMETIS_INC}"
                )
    else()
        message(STATUS "skipping parmetis, as zoltan2 doesn't need it.")
    endif()

    # check if we need to include metis
    check_defined(${Zoltan2_INCLUDE_DIRS}/Zoltan2_config.h HAVE_ZOLTAN2_METIS NEEDS_METIS)
    if(NOT ${NEEDS_METIS})
        message(STATUS "we need metis")
        target_include_directories(${Zoltan2LibraryTarget} SYSTEM INTERFACE
                "$ENV{METIS_INC}"
                )
    else()
        message(STATUS "skipping metis, as zoltan2 doesn't need it.")
    endif()

    message(STATUS "link with \${Zoltan2_LIBRARIES}: ${Zoltan2_LIBRARIES} ")
    message(STATUS "Zoltan2 support enabled")
else()
    message(STATUS "Zoltan2 support disabled")
endif()