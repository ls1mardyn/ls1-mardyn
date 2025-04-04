target_sources(MarDyn
    PRIVATE
        RealAccumVecTest.cpp
        RealCalcVecTest.cpp
        VCP1CLJRMMTest.cpp
    )

if(NOT ENABLE_AUTOPAS)
    target_sources(MarDyn PRIVATE VectorizedCellProcessorTest.cpp)
endif(NOT ENABLE_AUTOPAS)

