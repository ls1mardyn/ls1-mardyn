
$(info included tbbmalloc.mk!)

CXXFLAGS += -DTBB

ifneq ($(TBB_INC),)
  INCLUDES+= -I$(TBB_INC)
  $(info Appended TBB_INC)
else
  $(warning WARNING: TBB_INC not set)
endif

LDFLAGS += -ltbbmalloc_proxy -ltbbmalloc
ifneq ($(TBB_LIBDIR),)
  LDFLAGS += -L$(TBB_LIBDIR)
  $(info Appended TBB_LIBDIR)
else
  $(warning WARNING: TBB_LIBDIR not set)
endif

ifeq ($(CFG), cray-xt-gnu)
  LDFLAGS += -dynamic
  $(warning Using -dynamic for linking to tbbmallocproxy with cray-xt-gnu on Hazel Hen.)
endif