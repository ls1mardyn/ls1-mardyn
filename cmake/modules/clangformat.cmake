# Find clang-format
find_program(CLANG_FORMAT_BIN clang-format)

if (NOT CLANG_FORMAT_BIN)
    message(FATAL_ERROR "clang-format not found!")
endif()

# List of file extensions to format
set(FORMAT_EXTENSIONS cpp h hpp c)

# Recursively find all files with the specified extensions
file(GLOB_RECURSE ALL_SOURCE_FILES 
    ${PROJECT_SOURCE_DIR}/*.[ch]pp
    ${PROJECT_SOURCE_DIR}/*.[ch])

# Define clang-format target
add_custom_target(clang-format
    COMMAND ${CLANG_FORMAT_BIN} -i ${ALL_SOURCE_FILES}
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    COMMENT "Running clang-format"
)
