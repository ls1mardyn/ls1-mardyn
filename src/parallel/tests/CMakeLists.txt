target_sources(MarDyn
    PRIVATE
        DomainDecompBaseTest.cpp
        ZonalMethodTest.cpp
    )
if(ENABLE_MPI)
    target_sources(MarDyn
        PRIVATE
            CollectiveCommunicationTest.cpp
            CommunicationBufferTest.cpp
            DomainDecompositionTest.cpp
            KDDecompositionTest.cpp
            KDNodeTest.cpp
            NeighborAcquirerTest.cpp
            ZonalMethodTest.cpp
        )
endif(ENABLE_MPI)

