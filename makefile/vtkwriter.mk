
$(info included vtkwriter.mk!)

VTK_SOURCES = $(shell find ./ -name "*.cpp" | grep -v "/tests/" | grep "/vtk/")
VTK_INCDIR = /dss/dsshome1/lrz/sys/spack/release/19.2/opt/x86_avx2/xerces-c/3.2.1-gcc-jz3anmr/bin
VTK_LIBDIR = /dss/dsshome1/lrz/sys/spack/release/19.2/opt/x86_avx2/xerces-c/3.2.1-gcc-jz3anmr/lib

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
