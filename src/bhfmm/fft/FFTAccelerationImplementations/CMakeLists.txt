target_sources(MarDyn
    PRIVATE
        FFTAcceleration_2wayM2L.cpp
        FFTAcceleration_2wayM2L_fftw.cpp
        FFTAcceleration_2wayM2L_optFFT.cpp
        FFTAcceleration_matrices.cpp
        FFTAcceleration_matrices_fftw.cpp
        FFTAcceleration_matrices_optFFT.cpp
    )


add_subdirectory(block)

