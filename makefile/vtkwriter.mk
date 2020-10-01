
$(info included vtkwriter.mk!)

VTK_SOURCES = $(shell find ./ -name "*.cpp" | grep -v "/tests/" | grep "/vtk/")

SOURCES += $(VTK_SOURCES)
CXXFLAGS += -DVTK

INCLUDES += -isystem ../dependencies-external/libxsd
ifneq ($(VTK_INCDIR),)
  INCLUDES+= -isystem $(VTK_INCDIR)
  $(info Appended VTK_INCDIR)
else
  $(warning WARNING: VTK_INCDIR not set)
endif

LDFLAGS += -lxerces-c 
ifneq ($(VTK_LIBDIR),)
  LDFLAGS += -L$(VTK_LIBDIR)
  $(info Appended VTK_LIBDIR)
else
  $(warning WARNING: VTK_LIBDIR not set)
endif
