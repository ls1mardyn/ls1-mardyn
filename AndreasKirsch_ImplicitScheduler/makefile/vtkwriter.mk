
$(info included vtkwriter.mk!)

ifeq ($(VTK), 1)

VTK_SOURCES = $(shell find ./ -name "*.cpp" | grep -v "/tests/" | grep "/vtk/")

SOURCES += $(VTK_SOURCES)
INCLUDES += -I../dependencies-external/libxsd
LINKFLAGS += -lxerces-c
CXXFLAGS += -DVTK
endif