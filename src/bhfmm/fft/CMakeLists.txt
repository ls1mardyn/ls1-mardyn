target_sources(MarDyn
    PRIVATE
        FFTOrderReduction.cpp
        FFTSettings.cpp
    )


add_subdirectory(FFTAccelerationImplementations)
add_subdirectory(tools)
add_subdirectory(transferFunctionManager)

