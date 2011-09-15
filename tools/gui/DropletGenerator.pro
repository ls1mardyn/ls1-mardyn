TEMPLATE = lib

CONFIG += dll debug

MOC_DIR = moc_obj
OBJECTS_DIR = obj
DESTDIR=./libs

# ScenarioGenerator related headers
HEADERS  += src/Objects/Object.h
HEADERS  += src/Generators/Generator.h
HEADERS  += src/Parameters/Parameter.h
HEADERS  += src/Parameters/ParameterWithValue.h
HEADERS  += src/Parameters/ParameterWithIntValue.h
HEADERS  += src/Parameters/ParameterWithDoubleValue.h
HEADERS  += src/Parameters/ParameterWithChoice.h
HEADERS  += src/Parameters/ParameterCollection.h
HEADERS  += src/IO/WriteOutput.h
HEADERS  += src/Tokenize.h

# MD Generator headers
HEADERS  += src/Generators/md/MDGenerator.h
HEADERS  += src/Generators/md/PMFileReader.h
HEADERS  += src/Generators/md/ComponentParameters.h
HEADERS  += src/Generators/md/MardynConfiguration.h
HEADERS  += src/Generators/md/MardynConfigurationParameters.h
HEADERS  += src/Generators/md/MardynConfigLegacyWriter.h
HEADERS  += src/Generators/md/DropletGenerator.h
HEADERS  += src/Generators/md/DrawableMolecule.h
SOURCES  += src/Generators/md/MDGenerator.cpp
SOURCES  += src/Generators/md/PMFileReader.cpp
SOURCES  += src/Generators/md/ComponentParameters.cpp
SOURCES  += src/Generators/md/MardynConfiguration.cpp
SOURCES  += src/Generators/md/MardynConfigurationParameters.cpp
SOURCES  += src/Generators/md/DropletGenerator.cpp
SOURCES  += src/Generators/md/DrawableMolecule.cpp
SOURCES  += src/Generators/md/MardynConfigLegacyWriter.cpp

LIBS += -L. -lMardyn

INCLUDEPATH += ./src/
INCLUDEPATH += ../../src/
INCLUDEPATH += $(VTKINCLUDEPATH)