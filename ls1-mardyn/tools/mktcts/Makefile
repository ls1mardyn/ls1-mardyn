mkCP: main.o Domain.o Random.o
	g++ -g -o mkTcTS main.o Domain.o Random.o
main.o: main.cpp Domain.h
	g++ -Wall -g -c main.cpp
Domain.o: Domain.cpp Domain.h Random.h
	g++ -Wall -g -c Domain.cpp
Random.o: Random.cpp Random.h
	g++ -Wall -g -c Random.cpp

clean:
	rm -f *.o *~ mkTcTS *.bz2 *.xdr

tar: clean
	tar -cvjf mkTcTS.tar.bz2 *.h *.cpp Makefile

