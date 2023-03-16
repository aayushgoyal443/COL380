CC=g++
CFLAGS=-g -Wall
TARGET=main
TEST=0

compile_parallel:
	mpic++ -std=c++17 -o a2 main.cpp

run_parallel_parse:
	mpirun -np 4 ./a2 --taskid=1 --inputpath=test_cases/A2/test$(test)/test-input-$(test).gra --headerpath=test_cases/A2/test$(test)/test-header-$(test).dat --outputpath=results/result_$(test).txt  --verbose=0 --startk=$(start) --endk=$(end) --p=10

run_parallel_parse_verbose:
	mpirun -np 4 ./a2 --taskid=1 --inputpath=test_cases/A2/test$(test)/test-input-$(test).gra --headerpath=test_cases/A2/test$(test)/test-header-$(test).dat --outputpath=results/result_$(test)_verbose.txt  --verbose=1 --startk=$(start) --endk=$(end) --p=10