README

This code generates 10 million random elliptic Lambert's problems and uses them to compare the Exact solution to IvLam.

Step 1:
	- Compile the data generator and generate data.
		- Run these lines in your terminal:

gfortran -O3 -fopenmp generate_lambert.f90 -o generate_lambert
./generate_lambert

Step 2:
	- Compile the benchmark file and both functions, then run the test.
		- Run these lines in your terminal:


gfortran -O3 -march=native -ffree-line-length-none -fopenmp -cpp -DEXACT_STANDALONE ivLamRuntimeV2p50_739981p74401.f90 exactsolution.F90 benchmark_lambert.f90 -o bench_single.exe

.\bench_single.exe