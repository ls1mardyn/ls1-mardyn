$(info # Including cppunit.mk!)

# directory with the cppunit header files
#CPPUNIT_INC_DIR=/usr/include
#CXXFLAGS += -I$(CPPUNIT_INC_DIR)

# directory including the cppunit library
#CPPUNIT_LIB_DIR=/usr/lib64
#LDFLAGS += -L$(CPPUNIT_LIB_DIR)

#INCLUDES += -I../dependencies-external/cppunit-1.12.1/include

CPPUNIT_TESTS = $(shell find ./ -name "*.cpp" | grep -v "parallel/" | grep -v "vtk/" | grep "/tests/")
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

CXXFLAGS += -DUNIT_TESTS
LDFLAGS  += -ldl -lcppunit
SOURCES += $(CPPUNIT_TESTS)

