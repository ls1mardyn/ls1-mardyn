file(GLOB_RECURSE MY_SRC
        "*.cpp"
        "*.h"
        )

# if unit tests are disabled, remove the unit tests!
if(NOT ENABLE_UNIT_TESTS)
    list(FILTER MY_SRC EXCLUDE REGEX "/tests/")
endif()

# if mpi is not enabled, remove the uneeded source files
if(NOT ENABLE_MPI)
    # duplicate the list
    set(MY_SRC_BACK ${MY_SRC})
    # exclude everything from parallel
    list(FILTER MY_SRC EXCLUDE REGEX "/parallel/")

    # but include DomainDecompBase* and LoadCalc*
    list(FILTER MY_SRC_BACK INCLUDE REGEX "/parallel/")
    list(FILTER MY_SRC_BACK INCLUDE REGEX "DomainDecompBase|LoadCalc|Zonal")
    list(APPEND MY_SRC ${MY_SRC_BACK})
else()
    if(NOT ENABLE_ALLLBL)
        list(FILTER MY_SRC EXCLUDE REGEX "ALLLoadBalancer")
    endif()
endif()

# if fmm fft is disabled, remove everything in an fft directory
if(NOT ENABLE_FMM_FFT)
    list(FILTER MY_SRC EXCLUDE REGEX "/fft/")
endif()

# if vtk is disabled, remove everything in a vtk directory
if(NOT ENABLE_VTK)
    list(FILTER MY_SRC EXCLUDE REGEX "/vtk/")
endif()

if (NOT ENABLE_AUTOPAS)
    list(FILTER MY_SRC EXCLUDE REGEX "AutoPas")
else()
#    set(MY_SRC_BACK ${MY_SRC})
#    list(FILTER MY_SRC_BACK INCLUDE REGEX "AutoPas|ParticleIterator")
#    list(FILTER MY_SRC EXCLUDE REGEX "/particleContainer/")
#    list(APPEND MY_SRC ${MY_SRC_BACK})
    list(FILTER MY_SRC EXCLUDE REGEX "LinkedCells|VectorizedCellProcessorTest")
endif ()

# add the executable
ADD_EXECUTABLE(MarDyn
        ${MY_SRC})

# dependencies for lz4
if (ENABLE_LZ4)
    add_dependencies(MarDyn liblz4)
endif()

# we just add all libraries here. If a library is not set, it will simply be ignored.
TARGET_LINK_LIBRARIES(MarDyn
        ${BLAS_LIB}    # for armadillo
        ${LAPACK_LIB}  # for armadillo
        ${VTK_LIB}     # for VTK/xerces
        ${CPPUNIT_LIB} # for unit tests
        ${AUTOPAS_LIB} # for autopas
        ${LZ4_LIB}     # for LZ4 compression
        ${ALL_LIB}     # for ALL
        -Wl,--allow-multiple-definition # for MPI_TIMED linking
        )

ADD_TEST(
        NAME MarDyn_test
        COMMAND MarDyn -t -d ${PROJECT_SOURCE_DIR}/test_input
)
