compile:
	mpic++ -std=c++17 -o a2 main.cpp

run:
	mpirun -np 4 ./a2 test_cases/A2/test0/test-input-0.gra test_cases/A2/test0/test-header-0.dat test_cases/test0/output0.txt 1 10


compile_sample:
	mpic++ -std=c++17 -o sample sample.cpp
	mpirun -np 4 ./sample