target_sources(MarDyn
    PRIVATE
        Comp2Param.cpp
        Component.cpp
        FullMolecule.cpp
        MoleculeRMM.cpp
        Quaternion.cpp
    )

if(ENABLE_AUTOPAS)
    target_sources(MarDyn PRIVATE AutoPasSimpleMolecule.cpp)
endif(ENABLE_AUTOPAS)

add_subdirectory(mixingrules)

