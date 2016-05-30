PROJECT= mkSD
CC=g++
OBJECTS = main.o Component.o ConfigWriter.o GlobalStartGeometry.o  Molecule.o PhaseSpaceWriter.o RandomNumber.o
CFLAGS = -g -Wall -c 
LFLAGS = -g
#HEADERS = $(wildcard *.h)


$(PROJECT):$(OBJECTS)
	$(CC) -o $(PROJECT) $(OBJECTS) $(LFLAGS)

main.o: main.cpp Component.h Component.cpp  ConfigWriter.h ConfigWriter.cpp  PhaseSpaceWriter.h PhaseSpaceWriter.cpp 
	$(CC) $(CFLAGS) main.cpp 

Component.o: Component.h Component.cpp
	$(CC) $(CFLAGS)  Component.cpp

ConfigWriter.o: ConfigWriter.h ConfigWriter.cpp
	$(CC) $(CFLAGS) ConfigWriter.cpp

GlobalStartGeometry.o: GlobalStartGeometry.h GlobalStartGeometry.cpp RandomNumber.h RandomNumber.cpp
	$(CC) $(CFLAGS) GlobalStartGeometry.cpp 

Molecule.o: Molecule.h Molecule.cpp RandomNumber.h RandomNumber.cpp
	$(CC) $(CFLAGS) Molecule.cpp 

PhaseSpaceWriter.o: PhaseSpaceWriter.h PhaseSpaceWriter.cpp Component.h Component.cpp GlobalStartGeometry.h GlobalStartGeometry.cpp Molecule.h Molecule.cpp RandomNumber.h RandomNumber.cpp
	$(CC) $(CFLAGS) PhaseSpaceWriter.cpp 

RandomNumber.o: RandomNumber.h
	$(CC) $(CFLAGS) RandomNumber.cpp 
	
clean:
	rm -f *.o 
	rm -f mkSD

