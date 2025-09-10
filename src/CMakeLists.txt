# add the executable
if(MAMICO_COUPLING)
    add_library(MarDyn STATIC)
    target_include_directories(MarDyn SYSTEM PRIVATE "${MAMICO_SRC_DIR}")
    target_compile_definitions(MarDyn PUBLIC
        ${MAMICO_COMPILE_DEFINITIONS} ${MAMICO_MPI_DEFINITIONS}
        )
    if(MAMICO_ENABLE_FPIC)
        set_property(TARGET MarDyn PROPERTY POSITION_INDEPENDENT_CODE ON)
    endif()
else()
    add_executable(MarDyn)
endif()

target_sources(MarDyn
    PRIVATE
        Common.cpp
        Domain.cpp
        MarDyn.cpp
        Simulation.cpp
)

if(ENABLE_UNIT_TESTS)
    add_subdirectory(tests)
endif(ENABLE_UNIT_TESTS)


add_subdirectory(bhfmm)
add_subdirectory(ensemble)
add_subdirectory(integrators)
add_subdirectory(io)
add_subdirectory(longRange)
add_subdirectory(molecules)
add_subdirectory(particleContainer)
add_subdirectory(parallel)
add_subdirectory(plugins)
add_subdirectory(steereoCommands)
add_subdirectory(thermostats)
add_subdirectory(utils)


# we just add all libraries here. If a library is not set, it will simply be ignored.
target_link_libraries(MarDyn
        ${BLAS_LIB}    # for armadillo
        ${LAPACK_LIB}  # for armadillo
        ${VTK_LIB}     # for VTK/xerces
        ${CPPUNIT_LIB} # for unit tests
        ${AUTOPAS_LIB} # for autopas
        ${ADIOS2_LIB}  # for adios2
        ${LZ4_LIB}     # for LZ4 compression
        ${ALL_LIB}     # for ALL
        ${MPI_LIB}     # for MPI
        )

target_compile_definitions(MarDyn PUBLIC
        ${ADIOS2_COMPILE_DEFINITION}
)

add_test(
        NAME MarDyn_test
        COMMAND MarDyn -t -d ${PROJECT_SOURCE_DIR}/test_input
)
