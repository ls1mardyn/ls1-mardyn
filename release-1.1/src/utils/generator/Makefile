CC      = gcc
CXX     = g++
CFLAGS  = -g
CXFLAGS = $(CFLAGS)
LDFLAGS = -g

default: Generator_test Generator_test_c Lattice_test

%.o: %.cpp %.h
	$(CXX) $(CXXFLAGS) -c -o $@ $<


Lattice_test: Lattice_test.cpp Lattice.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) Lattice_test.cpp Lattice.o -o $@

Generator_test: Generator_test.o Generator.o Lattice.o Basis.o
	$(CXX) $^ $(LDFLAGS) -o $@

Generator_test_c: Generator_test_c.o Generator.o Lattice.o Basis.o
	$(CXX) $^ $(LDFLAGS) -o $@

.PHONY: clean

clean:
	$(RM) *.o
	$(RM) Generator_test Generator_test_c Lattice_test
