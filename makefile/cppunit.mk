
# this makefile builds the artefacts neccessary for mardyn to link against cppunit
# only variables used outside are CPPUNIT_INCLUDES and CPPUNIT_LINK_ARGS so far

#targets used outside are cppunit and cppunit_clean

$(warning included cppunit.mk!)

ifeq ($(TESTS), 1)

# sources of cppUnit
CPPUNIT_SOURCES = $(shell find ../dependencies-external/cppunit-1.12.1/src -name "*.cpp") $(shell find ../dependencies-external/cppunit-1.12.1/examples -name "*.cpp")

CPPUNIT_OBJECTS = $(CPPUNIT_SOURCES:.cpp=.o)

# if you want to use a cppunit-library which is already installed on your system,
# replace the following two lines to point to your installation and set CPPUNIT_OBJECTS empty
CPPUNIT_INCLUDES = -I../dependencies-external/cppunit-1.12.1/include
CPPUNIT_LINK_ARGS = $(CPPUNIT_OBJECTS) -ldl

# tests in "/parallel" are not consicered at the moment
CPPUNIT_TESTS = $(shell find ./ -name "*.cpp" | grep -v "parallel/" | grep "/tests/")
SOURCES += $(CPPUNIT_TESTS)
CXXFLAGS += -DUNIT_TESTS
CPPUNIT_TESTS_OBJECTS = $(CPPUNIT_TESTS:.cpp=.o)

cppunit: $(CPPUNIT_OBJECTS) $(CPPUNIT_TESTS_OBJECTS)

cppunit_clean:
	find ../dependencies-external/cppunit-1.12.1/src/ -type f -name '*.o' -delete
	find ../dependencies-external/cppunit-1.12.1/examples/ -type f -name '*.o' -delete

else
cppunit:
cppunit_clean:
endif
	
.PHONY: cppunit cppunit_clean
	
