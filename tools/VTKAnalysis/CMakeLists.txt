cmake_minimum_required(VERSION 3.3)

message(STATUS "cmake building vtk-analysis")

add_executable(vtk-analysis EXCLUDE_FROM_ALL
        main.cpp VTKFile.cpp VTKFile.h
        utils/xmlfileUnits.h
        utils/xmlfileUnits.cpp
        utils/xmlfile.h
        utils/xmlfile.cpp
        utils/OptionParser.h
        utils/OptionParser.cpp
        utils/FileUtils.h
        utils/FileUtils.cpp
        Molecule.h Analysis.h)
target_include_directories(vtk-analysis
        SYSTEM PUBLIC ../../libs/rapidxml
        )
TARGET_LINK_LIBRARIES(vtk-analysis
        ${MPI_LIB}     # for MPI
        )

