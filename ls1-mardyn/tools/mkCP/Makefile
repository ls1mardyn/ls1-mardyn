mkCP: main.o Graphit.o Domain.o Random.o
	g++ -g -o mkCP main.o Graphit.o Domain.o Random.o
main.o: main.cpp Domain.h Graphit.h
	g++ -Wall -g -c main.cpp
Graphit.o: Graphit.cpp Graphit.h Random.h
	g++ -Wall -g -c Graphit.cpp
Domain.o: Domain.cpp Domain.h Graphit.h Random.h
	g++ -Wall -g -c Domain.cpp
Random.o: Random.cpp Random.h
	g++ -Wall -g -c Random.cpp

clean:
	rm -f *.o *~ mkCP *.bz2 *.xdr

tar: clean
	tar -cvjf mkCP.tar.bz2 *.h *.cpp *.png *.txt Makefile
