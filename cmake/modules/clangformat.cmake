# Find clang-format
find_program(CLANG_FORMAT_BIN clang-format)

if (NOT CLANG_FORMAT_BIN)
    message(FATAL_ERROR "clang-format not found!")
endif()

# Recursively find all files with the specified extensions
file(GLOB_RECURSE ALL_SOURCE_FILES 
    ${PROJECT_SOURCE_DIR}/*.h
    ${PROJECT_SOURCE_DIR}/*.cpp
)

# Define clang-format target
add_custom_target(clang-format
    COMMAND ${CLANG_FORMAT_BIN} -i ${ALL_SOURCE_FILES}
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    COMMENT "Running clang-format"
)
