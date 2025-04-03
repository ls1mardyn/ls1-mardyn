target_sources(MarDyn
    PRIVATE
        AdaptivePseudoParticleContainer.cpp
        DttNode.cpp
        LeafNodesContainer.cpp
        ParticleCellPointers.cpp
        UniformPseudoParticleContainer.cpp
    )

if(ENABLE_UNIT_TESTS)
    add_subdirectory(tests)
endif(ENABLE_UNIT_TESTS)

