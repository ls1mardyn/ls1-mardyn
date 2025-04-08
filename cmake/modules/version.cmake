# this module runs git to determine the version and state of the code.
# it then edits MarDyn_version.h so that the version is available in C++

find_package(Git)
# we need git and the folder where git stores its info
if (Git_FOUND AND EXISTS "${MarDyn_SOURCE_DIR}/.git" )
    # branch name
    execute_process(
            COMMAND ${GIT_EXECUTABLE} rev-parse --abbrev-ref HEAD
            OUTPUT_VARIABLE MarDyn_VERSION_BRANCH)
    string(STRIP "${MarDyn_VERSION_BRANCH}" MarDyn_VERSION_BRANCH)
    # commit hash
    execute_process(
            COMMAND ${GIT_EXECUTABLE} rev-parse --short HEAD
            OUTPUT_VARIABLE MarDyn_VERSION_HASH)
    string(STRIP "${MarDyn_VERSION_HASH}" MarDyn_VERSION_HASH)
    # check dirty
    execute_process(
            COMMAND ${GIT_EXECUTABLE} diff --quiet --exit-code
            RESULT_VARIABLE is_dirty_result )
    if (${is_dirty_result})
        set(MarDyn_VERSION_IS_DIRTY "_dirty")
    else()
        set(MarDyn_VERSION_IS_DIRTY "")
    endif ()
else()
    message(WARNING "Could not find git or ${MarDyn_SOURCE_DIR}/.git folder! MARDYN_VERSION will be missing information.")
endif ()

configure_file(src/MarDyn_version.h.in src/MarDyn_version.h @ONLY)