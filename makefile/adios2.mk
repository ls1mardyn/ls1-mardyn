$(info included adios2.mk!)

ADIOS2_SOURCES = $(shell find ./ -name "*.cpp" | grep -v "/tests/" | grep "Adios")

SOURCES += $(ADIOS2_SOURCES)
CXXFLAGS += -DENABLE_ADIOS2

ifneq ($(ADIOS2_INCDIR),)
  INCLUDES+= -isystem $(ADIOS2_INCDIR)
  $(info Appended ADIOS2_INCDIR)
else
  INCLUDES+= $(shell adios2-config --cxx-flags)
endif
ifneq ($(ADIOS2_INCDIR),)
  $(error ADIOS2 was requested but no include paths could be found.)
endif

ifneq ($(ADIOS2_LIBDIR),)
  LDFLAGS += -L$(ADIOS2_LIBDIR)
  $(info Appended ADIOS2_LIBDIR)
else
  LDFLAGS += $(shell adios2-config --cxx-libs)
endif

ifneq ($(ADIOS2_LIBDIR),)
  $(error ADIOS2 was requested but no library paths could be found.)
endif
