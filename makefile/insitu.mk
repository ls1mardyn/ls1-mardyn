
$(info included insitu.mk!)

CXXFLAGS += -DINSITU

INCLUDES += -isystem ../dependencies-external/zeromq/include 
# ifneq ($(INSITU_INCDIR),)
#   INCLUDES+= -isystem $(INSITU_INCDIR)
#   $(info Appended INSITU_INCDIR)
# else
#   $(warning WARNING: INSITU_INCDIR not set)
# endif

LDFLAGS += -lzmq -L"/usr/lib64"
# ifneq ($(ZMQ_LIBDIR),)
#   LDFLAGS += -L$(ZMQ_LIBDIR)
#   $(info Appended ZMQ_LIBDIR)
# else
#   $(warning WARNING: ZMQ_LIBDIR not set)
# endif