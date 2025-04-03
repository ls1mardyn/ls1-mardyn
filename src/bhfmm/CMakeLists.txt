target_sources(MarDyn
    PRIVATE
        FastMultipoleMethod.cpp
    )

add_subdirectory(containers)
if(ENABLE_UNIT_TESTS)
    add_subdirectory(tests)
endif(ENABLE_UNIT_TESTS)


add_subdirectory(cellProcessors)
add_subdirectory(expansions)
add_subdirectory(pseudoParticles)
add_subdirectory(utils)

if(ENABLE_FMM_FFT)
    add_subdirectory(fft)
endif(ENABLE_FMM_FFT)

