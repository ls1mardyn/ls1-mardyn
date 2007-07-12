# machine type: SEQ or PAR
MACHTYPE = SEQ
# choose target: DEBUG|RELEASE
TARGET = RELEASE

SOURCES_SEQ = $(shell find ./ -name "*.cpp" | grep -v DomainDecomposition)
SOURCES_PAR = $(shell find ./ -name "*.cpp" | grep -v DomainDecompDummy)
SOURCES = $(SOURCES_$(MACHTYPE))

INCLUDES = -I$(shell pwd)

OBJECTS = $(shell echo $(SOURCES) | sed -e "s/\.cpp/.o/g")

BINARY = MarDyn

MAKEFILEDIR = $(dir $(realpath $(lastword $(MAKEFILE_LIST))))

# which compiler?
include $(MAKEFILEDIR)/cfg/gcc.in

CXX = $(CXX_$(MACHTYPE))
CXXFLAGS = $(CXXFLAGS_$(MACHTYPE)_$(TARGET))
LINKFLAGS = $(LINKFLAGS_$(MACHTYPE)_$(TARGET))


# only linker and domain decomposition need MPI
$(BINARY): $(OBJECTS)
	$(CXX) $(LINKFLAGS) -o $(BINARY) $(OBJECTS)

parallel/DomainDecomposition.o: parallel/DomainDecomposition.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# all other files can always be compiled with g++
%.o: %.cpp
	$(CXX_SEQ) $(CXXFLAGS) -c $< -o $@

.PHONY: help all clean cleanall depend

help:
	@echo "configuration: $(CFG_INFO)"
	@echo "targets:"
	@echo "make all	build $(BINARY)"
	@echo "make clean	delete object files"
	@echo "make cleanall	delete object files and executable"
	@echo "make depend	generate makefile dependency list"

all: $(BINARY)

clean:
	find -name '*.o' -exec rm -r {} \;

cleanall: clean
	rm -f $(BINARY)

depend:
	makedepend  -I$(INCLUDES) -- $(CXXFLAGS) -- $(SOURCES)

# DO NOT DELETE THIS LINE -- make depend depends on it.
