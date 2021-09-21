## This file tries to detect the vector with used by ARM SVE

set(SVE_VEC_WIDTH_OPTIONS "auto-detect;128;256;512")
set(SVE_VEC_WIDTH "auto-detect" CACHE STRING "Specifies the vector width (in bit) for SVE instructions (${SVE_VEC_WIDTH_OPTIONS}).")
# let ccmake and cmake-gui offer the options
set_property(CACHE SVE_VEC_WIDTH PROPERTY STRINGS ${SVE_VEC_WIDTH_OPTIONS})

if (SVE_VEC_WIDTH MATCHES "^auto-detect$")
    try_run(SVE_VEC_WIDTH_RUN_RES
            SVE_VEC_WIDTH_COMPILE_RES
            ${CMAKE_BINARY_DIR}
            ${MARDYN_SOURCE_DIR}/cmake/tests/sve_width_detection.cpp
            CMAKE_FLAGS "CMAKE_CXX_FLAGS=-march=native"
            COMPILE_OUTPUT_VARIABLE SVE_VEC_WIDTH_COMPILE_OUTPUT
            RUN_OUTPUT_VARIABLE SVE_VEC_WIDTH_RUN_OUTPUT
            )
    if(SVE_VEC_WIDTH_COMPILE_RES)
        if(SVE_VEC_WIDTH_RUN_RES MATCHES "^FAILED_TO_RUN$")
            message(FATAL_ERROR "SVE_VEC_WIDTH detector compiled successfully, but failed to run. Please set SVE_VEC_WIDTH manually!")
        else()
            message(STATUS "SVE_VEC_WIDTH detected to be ${SVE_VEC_WIDTH_RUN_OUTPUT}.")
            target_compile_definitions(MarDyn PUBLIC SVE_VEC_WIDTH=${SVE_VEC_WIDTH_RUN_OUTPUT})
        endif()
    else()
        message(STATUS "Could not determine SVE_VEC_WIDTH automatically. So not setting it.")
    endif()
else ()
    target_compile_definitions(MarDyn PUBLIC SVE_VEC_WIDTH=${SVE_VEC_WIDTH})
    message(STATUS "Setting SVE_VEC_WIDTH to ${SVE_VEC_WIDTH}.")
endif ()

