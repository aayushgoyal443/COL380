all: exec check

exec: main.o
	g++ -fopenmp -Wl,-rpath=. -g -o exec main.o -L/home/aayush/Desktop/COL380/Assignment1 -l:library.so

main.o: main.cpp library.hpp
	g++ -fopenmp -c main.cpp

check: check.o
	g++ -o check check.o

check.o: check.cpp
	g++ -c check.cpp

.PHONY: clean
clean:
	rm -f *.o exec check
