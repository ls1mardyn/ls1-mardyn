target_sources(MarDyn
    PRIVATE
        vtk-custom-decimal.cpp
        VTKGridCell.cpp
        VTKGridVertex.cpp
        VTKGridWriter.cpp
        VTKGridWriterImplementation.cpp
        VTKMoleculeWriter.cpp
        VTKMoleculeWriterImplementation.cpp
        vtk-punstructured.cpp
        vtk-unstructured.cpp
    )


if(ENABLE_UNIT_TESTS)
    add_subdirectory(tests)
endif(ENABLE_UNIT_TESTS)

