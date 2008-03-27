# machine type: SEQ or PAR
MACHTYPE = SEQ
# choose target: DEBUG|RELEASE
TARGET = RELEASE

SOURCES_SEQ = $(shell find ./ -name "*.cpp" | grep -v DomainDecomposition | grep -v External)
SOURCES_PAR = $(shell find ./ -name "*.cpp" | grep -v DomainDecompDummy | grep -v External)
SOURCES = $(SOURCES_$(MACHTYPE))

OBJECTS = $(shell echo $(SOURCES) | sed -e "s/\.cpp/.o/g")

# boost dependency
BOOST_VER = boost_1_34_0

LIBRARY_BOOST_LOCAL = $(shell find ../libs/$(BOOST_VER) -name 'libboost_program_options*.a' |grep release |grep threading-multi)
LIBRARY_BOOST_GLOBAL = $(shell find /usr/lib/ -name 'libboost_program_options*.a' |grep -e '-mt' ; find /usr/local/lib/ -name 'libboost_program_options*.a' |grep -e '-mt')


# tinyXPATH dependency
OBJECTS_TINYXPATH = $(patsubst %.cpp,%.o,$(shell find External/tinyxpath -name '*.cpp'))


INCLUDES = -I$(shell pwd) -I../libs/$(BOOST_VER)

BINARY = MarDyn

MAKEFILEDIR = $(dir $(realpath $(lastword $(MAKEFILE_LIST))))

# which compiler?
include $(MAKEFILEDIR)/cfg/gcc.in
#include ../makefile/cfg/gcc.in

CXX = $(CXX_$(MACHTYPE))
CXXFLAGS = $(CXXFLAGS_$(MACHTYPE)_$(TARGET))
LINKFLAGS = $(LINKFLAGS_$(MACHTYPE)_$(TARGET))


# only linker and domain decomposition need MPI
#$(BINARY): $(OBJECTS)
#	$(CXX) $(LINKFLAGS) -o $(BINARY) $(OBJECTS)
$(BINARY): $(OBJECTS) $(OBJECTS_TINYXPATH)
ifeq ($(strip $(LIBRARY_BOOST_LOCAL)),)
	$(CXX) $(LINKFLAGS) -o $(BINARY) $(OBJECTS) $(LIBRARY_BOOST_GLOBAL) $(OBJECTS_TINYXPATH)
else
	$(CXX) $(LINKFLAGS) -o $(BINARY) $(OBJECTS) $(LIBRARY_BOOST_LOCAL) $(OBJECTS_TINYXPATH)
endif

parallel/DomainDecomposition.o: parallel/DomainDecomposition.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# all other files can always be compiled with g++
%.o: %.cpp
	$(CXX_SEQ) $(CXXFLAGS) -c $< -o $@

.PHONY: help all clean cleanall depend check boost cleanboost

help:
	@echo "configuration: $(CFG_INFO)"
	@echo "targets:"
	@echo "make all		build $(BINARY)"
	@echo "make clean	delete object files"
	@echo "make cleanall	delete object files and executable"
	@echo "make depend	generate makefile dependency list"
	@echo "make check	check for needed libraries"
	@echo "make boost	build bundled boost library"
	@echo "make cleanboost	clean bundled boost library"

all: check $(BINARY)

clean:
	find -name '*.o' -exec rm -r {} \;

cleanall: clean
	rm -f $(BINARY)

depend:
	makedepend  -I$(INCLUDES) -- $(CXXFLAGS) -- $(SOURCES)

check:
ifeq ($(strip $(LIBRARY_BOOST_LOCAL)),)
ifeq ($(strip $(LIBRARY_BOOST_GLOBAL)),)
	@echo "ERROR: no boost library found. please run 'make boost' or install the boost 1.34.0 or newer package."
else
	@echo "SUCCESS: found libboost in $(shell dirname $(LIBRARY_BOOST_GLOBAL))"
endif
else
	@echo "SUCCESS: found libboost in $(shell dirname $(LIBRARY_BOOST_LOCAL))"
endif

boost:
	cd ../libs/$(BOOST_VER) && ./configure --with-libraries=program_options
	cd ../libs/$(BOOST_VER) && make

cleanboost:
	cd ../libs/$(BOOST_VER) && make distclean

# DO NOT DELETE THIS LINE -- make depend depends on it.
