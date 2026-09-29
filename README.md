# Exact Solution to Lambert's Problem: Benchmark vs. IvLam
 
This code generates **10 million random elliptic Lambert's problems** and uses them to compare the computational performance of the exact contour-integral solution to Lambert's problem [1, 2] against the IvLam solver [3, 4].
 
The exact solution recasts the root-finding problem for the universal variable *z* as a pole-finding problem. It then locates the pole with a ratio of complex contour integrals (via Cauchy's integral theorem), evaluated with the exponentially convergent composite trapezoidal rule, so no iteration is required. See [1] for the full derivation and [2] for the extended treatment, which also covers two-body dynamics, the Sundman transformation, and complex analysis.
 
IvLam is a complete Lambert solver based on the vercosine formulation. It uses an interpolated initial guess, built on custom KD-tree data structures, to seed a guarded root-solver that typically converges in 1–2 iterations, and it also provides exact first- and second-order sensitivities. See [3] for the method and [4] for the code.
 
## Requirements
 
- `gfortran` with OpenMP support
- Source files:
  - `generate_lambert.f90`: random Lambert problem generator
  - `exactsolution.F90`: exact contour-integral solver [1, 2]
  - `ivLamRuntimeV2p50_739981p74401.f90`: IvLam solver, v2.50 [3, 4]
  - `benchmark_lambert.f90`: benchmark driver
## Usage
 
### Step 1: Generate the test data
 
Compile the data generator and generate the data by running these lines in your terminal:
 
```bash
gfortran -O3 -fopenmp generate_lambert.f90 -o generate_lambert
./generate_lambert
```
 
### Step 2: Build and run the benchmark
 
Compile the benchmark file and both solvers, then run the test:
 
```bash
gfortran -O3 -march=native -ffree-line-length-none -fopenmp -cpp -DEXACT_STANDALONE ivLamRuntimeV2p50_739981p74401.f90 exactsolution.F90 benchmark_lambert.f90 -o bench_single.exe
```
 
On Windows:
 
```powershell
.\bench_single.exe
```
 
On Linux/macOS:
 
```bash
./bench_single.exe
```
 
## References
 
### Exact solution
 
1. Negrete, A., and Abdelkhalik, O. O., "Exact Solution to Lambert's Problem Using Contour Integrals," *Journal of Guidance, Control, and Dynamics*, Vol. 48, No. 4, 2025, pp. 771–780. https://doi.org/10.2514/1.G008499
2. Negrete Fernandez de Romarate, A., "An Exact Solution to Lambert's Problem," M.S. Thesis, Department of Aerospace Engineering, Iowa State University, Ames, IA, 2024.
### IvLam
 
3. Russell, R. P., "Complete Lambert Solver Including Second-Order Sensitivities," *Journal of Guidance, Control, and Dynamics*, Vol. 45, No. 2, 2022, pp. 196–212. https://doi.org/10.2514/1.G006089
4. Russell, R. P., "ivLam2 (v2.50)," Zenodo. https://zenodo.org/records/18102116
### BibTeX
 
```bibtex
@article{negrete2025exact,
  author  = {Negrete, Aimar and Abdelkhalik, Ossama O.},
  title   = {Exact Solution to Lambert's Problem Using Contour Integrals},
  journal = {Journal of Guidance, Control, and Dynamics},
  volume  = {48},
  number  = {4},
  pages   = {771--780},
  year    = {2025},
  doi     = {10.2514/1.G008499}
}
 
@mastersthesis{negrete2024thesis,
  author  = {Negrete Fernandez de Romarate, Aimar},
  title   = {An Exact Solution to Lambert's Problem},
  school  = {Iowa State University},
  address = {Ames, Iowa},
  year    = {2024}
}
 
@article{russell2022complete,
  author  = {Russell, Ryan P.},
  title   = {Complete Lambert Solver Including Second-Order Sensitivities},
  journal = {Journal of Guidance, Control, and Dynamics},
  volume  = {45},
  number  = {2},
  pages   = {196--212},
  year    = {2022},
  doi     = {10.2514/1.G006089}
}
 
@misc{russell_ivlam2_v250,
  author       = {Russell, Ryan P.},
  title        = {ivLam2 (v2.50)},
  howpublished = {Zenodo},
  url          = {https://zenodo.org/records/18102116}
}
```
