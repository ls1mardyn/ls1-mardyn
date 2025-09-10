target_sources(MarDyn
    PRIVATE
        BoxDomain.cpp
        CanonicalEnsemble.cpp
        CavityEnsemble.cpp
        ChemicalPotential.cpp
        EnsembleBase.cpp
        GrandCanonicalEnsemble.cpp
        PressureGradient.cpp
    )

if(ENABLE_UNIT_TESTS)
    add_subdirectory(tests)
endif(ENABLE_UNIT_TESTS)

