target_sources(MarDyn
    PRIVATE
        ASCIIReader.cpp
        BinaryReader.cpp
        CavityWriter.cpp
        CheckpointWriter.cpp
        CommunicationPartnerWriter.cpp
        CubicGridGeneratorInternal.cpp
        DecompWriter.cpp
        EnergyLogWriter.cpp
        FlopRateWriter.cpp
        GammaWriter.cpp
        HaloParticleWriter.cpp
        IOHelpers.cpp
        KDTreePrinter.cpp
        LoadBalanceWriter.cpp
        MaxWriter.cpp
        MemoryProfiler.cpp
        Mkesfera.cpp
        MmpldWriter.cpp
        MmspdBinWriter.cpp
        MmspdWriter.cpp
        MPICheckpointWriter.cpp
        MultiObjectGenerator.cpp
        ObjectGenerator.cpp
        ODF.cpp
        PerCellGenerator.cpp
        PovWriter.cpp
        RDF.cpp
        ReplicaGenerator.cpp
        ResultWriter.cpp
        SysMonOutput.cpp
        TaskTimingProfiler.cpp
        TcTS.cpp
        TimerProfiler.cpp
        TimerWriter.cpp
        XyzWriter.cpp
    )

if(ENABLE_ADIOS2)
    target_sources(MarDyn PRIVATE Adios2Reader.cpp Adios2Writer.cpp)
endif(ENABLE_ADIOS2)

if(ENABLE_UNIT_TESTS)
    add_subdirectory(tests)
endif(ENABLE_UNIT_TESTS)


if(ENABLE_VTK)
    add_subdirectory(vtk)
endif(ENABLE_VTK)

