target_sources(MarDyn
    PRIVATE
        Comp2Param.cpp
        Component.cpp
        FullMolecule.cpp
        MoleculeInterface.cpp
        MoleculeRMM.cpp
        Quaternion.cpp
    )

if(ENABLE_AUTOPAS)
    target_sources(MarDyn PRIVATE AutoPasSimpleMolecule.cpp)
endif(ENABLE_AUTOPAS)

if(ENABLE_UNIT_TESTS)
    add_subdirectory(tests)
endif(ENABLE_UNIT_TESTS)


add_subdirectory(mixingrules)

