
$(info included vtkwriter.mk!)

VTK_SOURCES = $(shell find ./ -name "*.cpp" | grep -v "/tests/" | grep "/vtk/")

SOURCES += $(VTK_SOURCES)
INCLUDES += -I../dependencies-external/libxsd
LDFLAGS += -lxerces-c
CXXFLAGS += -DVTK
