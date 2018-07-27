MACHTYPE = SEQ
TARGET = RELEASE

INCLUDES = -I$(shell pwd)

include ../makefile/cfg/gcc.in

CXX = $(CXX_$(MACHTYPE))
CXXFLAGS = $(CXXFLAGS_$(MACHTYPE)_$(TARGET))
LINKFLAGS = $(LINKFLAGS_$(MACHTYPE)_$(TARGET))

OBJECTS_TINYXPATH = $(patsubst %.cpp,%.o,$(shell find ../src/External/tinyxpath -name '*.cpp'))

all: moldy2mardyn mdproject2mardyn

moldy2mardyn: moldy2mardyn.o $(OBJECTS_TINYXPATH)
	$(CXX) $(LINKFLAGS) -o moldy2mardyn moldy2mardyn.o $(OBJECTS_TINYXPATH)
mdproject2mardyn: mdproject2mardyn.o $(OBJECTS_TINYXPATH)
	$(CXX) $(LINKFLAGS) -o mdproject2mardyn mdproject2mardyn.o $(OBJECTS_TINYXPATH)

%.o: %.cpp
	$(CXX_SEQ) $(CXXFLAGS) -c $< -o $@

clean:
	find -name '*.o' -exec rm -f {} \;
	rm -f mdproject2mardyn moldy2mardyn
