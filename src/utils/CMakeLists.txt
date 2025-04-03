target_sources(MarDyn
    PRIVATE
        Coordinate3D.cpp
        Expression.cpp
        FileUtils.cpp
        Logger.cpp
        Math.cpp
        MPI_Info_object.cpp
        nnls.cpp
        OptionParser.cpp
        PrintThreadPinningToCPU.cpp
        Random.cpp
        Region.cpp
        SteereoIntegration.cpp
        SysMon.cpp
        Testing.cpp
        TestWithSimulationSetup.cpp
        Timer.cpp
        xmlfile.cpp
        xmlfileUnits.cpp
    )

if(ENABLE_UNIT_TESTS)
    add_subdirectory(tests)
endif(ENABLE_UNIT_TESTS)


add_subdirectory(generator)

