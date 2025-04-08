target_sources(MarDyn
    PRIVATE
        fft_utils.cpp
        FFTW_Helper.cpp
    )


add_subdirectory(optimizedFFT)

