$(info included adios2.mk!)

ADIOS2_SOURCES = $(shell find ./ -name "*.cpp" | grep -v "/tests/" | grep "Adios2")

SOURCES += $(ADIOS2_SOURCES)
CXXFLAGS += -DENABLE_ADIOS2

ifneq ($(ADIOS2_CONFIG_BIN),)
  ADIOS2_CONFIG_BIN="adios2-config"
endif

CXXFLAGS += $(shell $(ADIOS2_CONFIG_BIN) --cxx-flags)
LDFLAGS += $(shell $(ADIOS2_CONFIG_BIN) --cxx-libs)