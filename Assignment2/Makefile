CC=g++
CFLAGS=-g -Wall
TARGET=main
TEST=0

compile_parallel:
	mpic++ -std=c++17 -o a2 main.cpp

run_parallel_parse:
	mpirun -np 4 ./a2 --taskid=1 --inputpath=test_cases/A2/test0/test-input-0.gra --headerpath=test_cases/A2/test0/test-header-0.dat --outputpath=test_cases/A2/test0/output0.txt  --verbose=1 --startk=1 --endk=10 --p=10

run_parallel:
	mpirun -np 4 ./a2 test_cases/A2/test$(test)/test-input-$(test).gra test_cases/A2/test$(test)/test-header-$(test).dat test_cases/A2/test$(test)/output$(test).txt $(start) $(end)

compile_serial:
	g++ -std=c++17 -o a2_serial main_serial.cpp

run_serial:
	./a2_serial our_ass/test_cases/test$(test)/test-input-$(test).gra our_ass/test_cases/test$(test)/test-header-$(test).dat our_ass/test_cases/test$(test)/output$(test).txt $(start) $(end)
