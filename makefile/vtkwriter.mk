
$(info included vtkwriter.mk!)

VTK_SOURCES = $(shell find ./ -name "*.cpp" | grep -v "/tests/" | grep "/vtk/")

SOURCES += $(VTK_SOURCES)
CXXFLAGS += -DVTK
INCLUDES += -I../dependencies-external/libxsd -I$(VTK_INCDIR)
LDFLAGS += -lxerces-c -L$(VTK_LIBDIR)
