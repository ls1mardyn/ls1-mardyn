$(info # Including cppunit.mk!)

# directory with the cppunit header files
#CPPUNIT_INC_DIR=/usr/include
#CXXFLAGS += -I$(CPPUNIT_INC_DIR)

# directory including the cppunit library
#CPPUNIT_LIB_DIR=/usr/lib64
#LDFLAGS += -L$(CPPUNIT_LIB_DIR)


# tests in "/parallel" are not consicered at the moment
CPPUNIT_TESTS = $(shell find ./ -name "*.cpp" | grep -v "parallel/" | grep -v "vtk/" | grep "/tests/")
ifeq ($(VTK), 1)
CPPUNIT_TESTS += $(shell find ./ -name "*.cpp" | grep -v "parallel/" | grep "vtk/tests/")
endif

CXXFLAGS += -DUNIT_TESTS
LDFLAGS  += -ldl -lcppunit
SOURCES += $(CPPUNIT_TESTS)


.PHONY: test
	
ifeq ($(PARTYPE), PAR)
$(info # Parallel defined!) 
MPICMD=mpirun -n 2 
endif

test: $(BINARY)
	$(info # Running test with $(MPICMD)) 
	$(MPICMD) ./$(BINARY) -t -d ../test_input/

