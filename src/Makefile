# machine type: SEQ or PAR
MACHTYPE = PAR
# choose target: DEBUG|RELEASE
#TARGET = RELEASE
TARGET = DEBUG
OBJECTS_PAR = $(shell find ./ -name "*.cpp" | sed -e "s/\.cpp/.o/g" | grep -v DomainDecompDummy)
OBJECTS_SEQ = $(shell find ./ -name "*.cpp" | sed -e "s/\.cpp/.o/g" | grep -v 'DomainDecomposition\|KDDecomposition')

OBJECTS = $(OBJECTS_$(MACHTYPE))

BINARY = MarDyn

CXX_SEQ = g++
CXX_PAR = mpicxx
CXXFLAGS_SEQ_DEBUG = -I$(shell pwd) -Wall -DDEBUG -g 
CXXFLAGS_SEQ_RELEASE = -I$(shell pwd)  -O2 -DNDEBUG 
LINKFLAGS_SEQ_DEBUG = -O0 -g
LINKFLAGS_SEQ_RELEASE = 

CXXFLAGS_PAR_DEBUG = -I$(shell pwd) -Wall -g -DDEBUG -DPARALLEL
CXXFLAGS_PAR_RELEASE = -I$(shell pwd) -O2 -DNDEBUG -DPARALLEL
LINKFLAGS_PAR_DEBUG = -O0 -g -DPARALLEL
LINKFLAGS_PAR_RELEASE = -DPARALLEL


CXX = $(CXX_$(MACHTYPE))
CXXFLAGS = $(CXXFLAGS_$(MACHTYPE)_$(TARGET))
LINKFLAGS = $(LINKFLAGS_$(MACHTYPE)_$(TARGET)) -lrt


# only linker and domain decomposition need MPI
$(BINARY): $(OBJECTS)
	$(CXX) $(LINKFLAGS) -o $(BINARY) $(OBJECTS)

parallel/DomainDecomposition.o: parallel/DomainDecomposition.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@
parallel/KDDecomposition.o: parallel/KDDecomposition.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@
Simulation.o: Simulation.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# all other files can always be compiled with g++
%.o: %.cpp
	$(CXX_SEQ) $(CXXFLAGS) -c $< -o $@

clean:
	find -name '*.o' -exec rm -r {} \;

cleanall: clean
	rm -f $(BINARY)
