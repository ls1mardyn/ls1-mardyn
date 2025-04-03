target_sources(MarDyn
    PRIVATE
        DomainDecompBase.cpp
        ForceHelper.cpp
        LoadCalc.cpp
    )

if(ENABLE_MPI)
    target_sources(MarDyn
        PRIVATE
            CommunicationBuffer.cpp
            CommunicationPartner.cpp
            DomainDecompMPIBase.cpp
            DomainDecomposition.cpp
            GeneralDomainDecomposition.cpp
            KDDecomposition.cpp
            KDDStaticValues.cpp
            KDNode.cpp
            MPIKDNode.cpp
            NeighborAcquirer.cpp
            NeighbourCommunicationScheme.cpp
            NonBlockingMPIHandlerBase.cpp
            NonBlockingMPIMultiStepHandler.cpp
            ParticleDataFull.cpp
            ParticleDataRMM.cpp
            ParticleForceData.cpp
            ResilienceComm.cpp
            StaticIrregDomainDecomposition.cpp
        )
    if(ENABLE_ALLLBL)
        target_sources(MarDyn PRIVATE ALLLoadBalancer.cpp)
    endif(ENABLE_ALLLBL)
endif(ENABLE_MPI)

if(ENABLE_UNIT_TESTS)
    add_subdirectory(tests)
endif(ENABLE_UNIT_TESTS)


add_subdirectory(boundaries)
add_subdirectory(ZonalMethods)

