$(info # Including cppunit.mk!)

#CPPUNIT_DIR=../dependencies-external/cppunit-1.12.1

CXXFLAGS += -DUNIT_TESTS

ifeq ($(CPPUNIT_DIR),)
ifneq ($(CPPUNIT_BASE),)
CPPUNIT_DIR=$(CPPUNIT_BASE)
endif
endif

ifneq ($(CPPUNIT_DIR),)
CPPUNIT_LIB_DIR=$(CPPUNIT_DIR)/lib
CPPUNIT_INC_DIR=$(CPPUNIT_DIR)/include
CXXFLAGS += -isystem $(CPPUNIT_INC_DIR)
LDFLAGS += -L$(CPPUNIT_LIB_DIR)
endif

LDFLAGS  += -lcppunit


CPPUNIT_TESTS = $(shell find ./ -name "*.cpp" | grep -v "AutoPas" | grep -v "parallel/" | grep -v "vtk/" | grep "/tests/")
ifneq ($(PARTYPE), PAR)
#include the sequential DomainDecompBaseTest (if SEQTYPE == PAR, it will be included below with the parallel tests)
CPPUNIT_TESTS += $(shell find ./ -name "*.cpp" | grep "parallel/tests/DomainDecompBaseTest")
endif
ifeq ($(VTK), 1)
CPPUNIT_TESTS += $(shell find ./ -name "*.cpp" | grep -v "parallel/" | grep "vtk/tests/")
endif
ifeq ($(PARTYPE), PAR)
$(info ADDING PARALLEL TESTS!)
CPPUNIT_TESTS += $(shell find ./ -name "*.cpp" | grep "parallel/tests/")
endif

SOURCES += $(CPPUNIT_TESTS)

