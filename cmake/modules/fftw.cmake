include(FetchContent)

# Try to find FFTW in the system
find_library(FFTW_LIB fftw3)

# If not found, download and build from GitHub
if(NOT FFTW_LIB)
    message(STATUS "FFTW not found. Downloading and building from source...")

    FetchContent_Declare(
            fftw
            GIT_REPOSITORY https://github.com/FFTW/fftw3.git
            GIT_TAG        master  # Change to stable tag if needed, e.g., fftw-3.3.10
    )

    # Fetch and build FFTW
    FetchContent_MakeAvailable(fftw)

    # Manually set FFTW_LIB if needed
    if(TARGET fftw)
        set(FFTW_LIB fftw)
    else()
        message(FATAL_ERROR "Failed to build FFTW from source!")
    endif()
else()
    message(STATUS "Found FFTW: ${FFTW_LIB}")
endif()

# Provide FFTW_LIB to the main CMakeLists
set(FFTW_LIB ${FFTW_LIB} PARENT_SCOPE)
