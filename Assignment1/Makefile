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

run1:
	make clean && make && ./exec input1 out1

run_check1:
	make clean && make && ./exec input1 out1 && ./check out1 output1

run2:
	make clean && make && ./exec input2 out2

run_check2:
	make clean && make && ./exec input2 out2 && ./check out2 output2