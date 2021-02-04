TEMPLATE = lib

#CONFIG += dll debug
CONFIG += staticlib debug mardyn
lessThan(QT_MAJOR_VERSION, 5) {
	QMAKE_CXXFLAGS += -std=c++17
} else {
	CONFIG += c++17
}
QMAKE_CXXFLAGS += -O3
MOC_DIR = moc_obj
OBJECTS_DIR = obj

mardyn {
DESTDIR=./staticlibs
} else {
DESTDIR=./libs
}

# ScenarioGenerator related headers
HEADERS  += src/Generators/Generator.h
HEADERS  += src/Parameters/Parameter.h
HEADERS  += src/Parameters/ParameterWithValue.h
HEADERS  += src/Parameters/ParameterWithIntValue.h
HEADERS  += src/Parameters/ParameterWithLongIntValue.h
HEADERS  += src/Parameters/ParameterWithDoubleValue.h
HEADERS  += src/Parameters/ParameterWithChoice.h
HEADERS  += src/Parameters/ParameterCollection.h
HEADERS  += src/IO/WriteOutput.h
HEADERS  += src/Tokenize.h

mardyn {
SOURCES  += src/Parameters/Parameter.cpp
SOURCES  += src/Parameters/ParameterWithValue.cpp
SOURCES  += src/Parameters/ParameterWithIntValue.cpp
SOURCES  += src/Parameters/ParameterWithLongIntValue.cpp
SOURCES  += src/Parameters/ParameterWithDoubleValue.cpp
SOURCES  += src/Parameters/ParameterWithBool.cpp
SOURCES  += src/Parameters/ParameterWithStringValue.cpp
SOURCES  += src/Parameters/ParameterWithChoice.cpp
SOURCES  += src/Parameters/ParameterCollection.cpp
SOURCES  += src/IO/WriteOutput.cpp
SOURCES  += src/Tokenize.cpp

DEFINES += MARDYN
} else {
HEADERS  += src/Objects/Object.h
greaterThan(QT_MAJOR_VERSION, 4): QT += widgets core gui
INCLUDEPATH += $(VTKINCLUDEPATH)
LIBS += -L. -lMardyn
HEADERS  += generators/common/DrawableMolecule.h
SOURCES  += generators/common/DrawableMolecule.cpp
}

# MD Generator headers
HEADERS  += generators/MDGenerator.h
HEADERS  += generators/AqueousNaClGenerator.h
HEADERS  += generators/common/ComponentParameters.h
HEADERS  += generators/common/MardynConfiguration.h
HEADERS  += generators/common/MardynConfigurationParameters.h
HEADERS  += generators/common/PMFileReader.h
HEADERS  += generators/common/MardynConfigLegacyWriter.h
HEADERS  += generators/common/eig3.h
HEADERS  += generators/common/PrincipalAxisTransform.h
HEADERS  += generators/common/OutputConfiguration.h

SOURCES  += generators/MDGenerator.cpp
SOURCES  += generators/AqueousNaClGenerator.cpp
SOURCES  += generators/common/ComponentParameters.cpp
SOURCES  += generators/common/MardynConfiguration.cpp
SOURCES  += generators/common/MardynConfigurationParameters.cpp
SOURCES  += generators/common/PMFileReader.cpp
SOURCES  += generators/common/MardynConfigLegacyWriter.cpp
SOURCES  += generators/common/eig3.cpp
SOURCES  += generators/common/PrincipalAxisTransform.cpp
SOURCES  += generators/common/OutputConfiguration.cpp

INCLUDEPATH += ./src/
INCLUDEPATH += ../../src/
INCLUDEPATH += ../../src/External/
INCLUDEPATH += ../../libs/armadillo/
INCLUDEPATH += ../../libs/rapidxml/
