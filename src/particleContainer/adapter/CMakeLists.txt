target_sources(MarDyn
    PRIVATE
        FlopCounter.cpp
        LegacyCellProcessor.cpp
        ODFCellProcessor.cpp
        RDFCellProcessor.cpp
        VCP1CLJRMM.cpp
        VectorizedCellProcessor.cpp
    )

if(ENABLE_UNIT_TESTS)
    add_subdirectory(tests)
endif(ENABLE_UNIT_TESTS)

