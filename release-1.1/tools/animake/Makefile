mkCP: main.o Domain.o Random.o
	g++ -g -o ./animake main.o Domain.o Random.o
main.o: main.cpp Domain.h
	g++ -Wall -g -c main.cpp
Domain.o: Domain.cpp Domain.h Random.h
	g++ -Wall -g -c Domain.cpp
Random.o: Random.cpp Random.h
	g++ -Wall -g -c Random.cpp

clean:
	rm -f *.o *~ animake *.bz2 *.xdr

tar: clean
	tar -cvjf animake.tar.bz2 *.h *.cpp *.txt *.ods Makefile

