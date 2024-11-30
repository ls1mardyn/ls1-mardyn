grideval: grideval.o Domain.o
	g++ -g -o ./grideval grideval.o Domain.o
grideval.o: grideval.cpp grideval.h Domain.h
	g++ -Wall -g -c grideval.cpp
Domain.o: Domain.cpp Domain.h
	g++ -Wall -g -c Domain.cpp

clean:
	rm -f *.o *~ grideval *.bz2 *.xdr

tar: clean
	tar -cvjf grideval.tar.bz2 *.h *.cpp Makefile
