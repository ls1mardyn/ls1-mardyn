
$(info included fftw.mk!)

# TODO: this includes all fft sources, not just only the fftw-requiring ones
FFTW_SOURCES = $(shell find ./ -name "*.cpp" | grep "/fft/")

SOURCES += $(FFTW_SOURCES)
CXXFLAGS += -DFFTW

ifneq ($(FFTW_INCDIR),)
  INCLUDES+= -I$(FFTW_INCDIR)
  $(info Appended FFTW_INCDIR)
else
  $(warning WARNING: FFTW_INCDIR not set)
endif

LDFLAGS += -lfftw3 -lfftw3f
ifneq ($(FFTW_LIBDIR),)
  LDFLAGS += -L$(FFTW_LIBDIR)
  $(info Appended FFTW_LIBDIR)
else
  $(warning WARNING: FFTW_LIBDIR not set)
endif
