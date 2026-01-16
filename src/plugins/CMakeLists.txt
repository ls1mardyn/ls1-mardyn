target_sources(MarDyn
    PRIVATE
        COMaligner.cpp
        DirectedPM.cpp
        Dropaccelerator.cpp
        Dropaligner.cpp
        EnergyRAPL.cpp
        ExamplePlugin.cpp
        FixRegion.cpp
        PinningInfo.cpp
        InMemoryCheckpointing.cpp
        LoadImbalanceThroughSleepPlugin.cpp
        MamicoCoupling.cpp
        MaxCheck.cpp
        Mirror.cpp
        MirrorSystem.cpp
        Permittivity.cpp
        PluginFactory.cpp
        SpatialProfile.cpp
        VectorizationTuner.cpp
        WallPotential.cpp
    )

if(ENABLE_UNIT_TESTS)
    add_subdirectory(tests)
endif(ENABLE_UNIT_TESTS)


add_subdirectory(NEMD)
add_subdirectory(profiles)

