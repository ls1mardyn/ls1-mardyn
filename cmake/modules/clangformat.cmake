
set(INCLUDE_DIRS
        "${PROJECT_SOURCE_DIR}/src/"
)

# Find clang-format
find_program(CLANG_FORMAT NAMES clang-format)

if(CLANG_FORMAT)
    execute_process(
            COMMAND ${CLANG_FORMAT} --version
            OUTPUT_VARIABLE clang_format_version
            OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    if(clang_format_version MATCHES "version 18\\.[0-9]+\\.[0-9]+")
        message(STATUS "Found clang-format version 18: ${CLANG_FORMAT}")

        file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/.clang-format-stamps")

        # Create a custom target for formatting
#        add_custom_target(clangformat
#                COMMAND ${CMAKE_COMMAND} -E echo "Running clang-format..."
#        )

        # Find all source files
        set(ALL_SOURCE_FILES)
        foreach(DIR ${INCLUDE_DIRS})
            file(GLOB_RECURSE DIR_SOURCE_FILES
                    "${DIR}/*.cpp"
                    "${DIR}/*.h"
            )
            list(APPEND ALL_SOURCE_FILES ${DIR_SOURCE_FILES})
        endforeach()

        # Create individual file targets
        set(STAMP_FILES)

        foreach(_file ${ALL_SOURCE_FILES})
            # Generate unique stamp file name
            string(MD5 FILE_HASH "${_file}")
            get_filename_component(FILE_NAME "${_file}" NAME)
            set(STAMP_FILE "${CMAKE_BINARY_DIR}/.clang-format-stamps/${FILE_NAME}_${FILE_HASH}.stamp")

            add_custom_command(
                    OUTPUT ${STAMP_FILE}
                    COMMAND ${CLANG_FORMAT} -style=file -i "${_file}"
                    COMMAND ${CMAKE_COMMAND} -E touch "${STAMP_FILE}"
                    DEPENDS "${_file}"
                    COMMENT "Formatting ${_file}"
            )

            list(APPEND STAMP_FILES ${STAMP_FILE})
        endforeach()

#            if(SOURCE_FILES)
#                add_custom_command(
#                        TARGET clangformat
#                        COMMAND ${CLANG_FORMAT} -style=file -i ${SOURCE_FILES}
#                        COMMENT "Formatting files in ${DIR}"
#                )
#            else()
#                message(STATUS "No source files found in ${DIR}")
#            endif()
#
#            unset(SOURCE_FILES)  # Clear for next iteration
#        endforeach()
        # Create main formatting target
        add_custom_target(clangformat
                DEPENDS ${STAMP_FILES}
                COMMENT "Formatting changed files..."
        )
        message(STATUS "clang-format target added")
    else()
        message(WARNING "clang-format-18 not found. Found version: ${clang_format_version}")
        set(CLANG_FORMAT CLANG_FORMAT-NOTFOUND)
    endif()
else()
    message(STATUS "clang-format-18 not found")
endif()


## Recursively find all files with the specified extensions
#file(GLOB_RECURSE ALL_SOURCE_FILES
#    ${PROJECT_SOURCE_DIR}/*.h
#    ${PROJECT_SOURCE_DIR}/*.cpp
#)
#
## Define clang-format target
#add_custom_target(clang-format
#    COMMAND ${CLANG_FORMAT_BIN} -i ${ALL_SOURCE_FILES}
#    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
#    COMMENT "Running clang-format"
#)
