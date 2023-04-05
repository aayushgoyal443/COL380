CC=g++
CFLAGS=-g -Wall
TARGET=main
TEST=0

compile:
	mpic++ -O2 -fopenmp -w -std=c++17 -o a3 main_dsu.cpp

run_task1:
	mpirun -np 4 ./a3 --taskid=1 --inputpath=test_cases/A3/test$(test)/test-input-$(test).gra --headerpath=test_cases/A3/test$(test)/test-header-$(test).dat --outputpath=results/task1/result_$(test).txt --verbose=0 --startk=$(start) --endk=$(end) --p=$(p)

run_task1_verbose:
	mpirun -np 4 ./a3 --taskid=1 --inputpath=test_cases/A3/test$(test)/test-input-$(test).gra --headerpath=test_cases/A3/test$(test)/test-header-$(test).dat --outputpath=results/task1/result_$(test)_verbose.txt --verbose=1 --startk=$(start) --endk=$(end) --p=$(p)


run_task2:
	mpirun -np 4 ./a3 --taskid=2 --inputpath=test_cases/A3/test$(test)/test-input-$(test).gra --headerpath=test_cases/A3/test$(test)/test-header-$(test).dat --outputpath=results/task2/result_$(test).txt  --verbose=0 --startk=$(start) --endk=$(end) --p=$(p)

run_task2_verbose:
	mpirun -np 4 ./a3 --taskid=2 --inputpath=test_cases/A3/test$(test)/test-input-$(test).gra --headerpath=test_cases/A3/test$(test)/test-header-$(test).dat --outputpath=results/task2/result_$(test)_verbose.txt  --verbose=1 --startk=$(start) --endk=$(end) --p=$(p)

clean: 
	rm a3 results/task1/* results/task2/*