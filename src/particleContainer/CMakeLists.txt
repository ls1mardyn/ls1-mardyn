target_sources(MarDyn
    PRIVATE
        FullParticleCell.cpp
        ParticleCellBase.cpp
        ParticleCellRMM.cpp
        ParticleContainer.cpp
    )

if(ENABLE_AUTOPAS)
    target_sources(MarDyn PRIVATE AutoPasContainer.cpp)
else()
    target_sources(MarDyn PRIVATE LinkedCells.cpp)
endif(ENABLE_AUTOPAS)

if(ENABLE_UNIT_TESTS)
    add_subdirectory(tests)
endif(ENABLE_UNIT_TESTS)


add_subdirectory(adapter)

if(ENABLE_AUTOPAS)
    add_subdirectory(AutoPasTemplateInstantiations)
endif(ENABLE_AUTOPAS)

