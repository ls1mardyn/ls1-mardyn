target_sources(MarDyn
    PRIVATE
        ParticleCellRMMTest.cpp
        ParticleContainerFactory.cpp
        ParticleContainerTest.cpp
    )

if(ENABLE_AUTOPAS)
    target_sources(MarDyn PRIVATE AutoPasContainerTest.cpp)
else()
    target_sources(MarDyn PRIVATE LinkedCellsTest.cpp)
endif(ENABLE_AUTOPAS)

