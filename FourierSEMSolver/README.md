# Fourier-SEM solver
All code written by Nikolas Borrel-Jensen (2023).

Code for reproducing results from 'Accelerated sound propagation simulations using an error-free Fourier method coupled with the spectral element method' (N. Borrel-Jensen, A. P. Engsig-Karup, M. Hornikx, C.H. Jeong, 2022).

* `mainCoupledMethods.m`: Main functionality for generating data and plots for a variety of couplings (`.FDTD_FDTD` | `.FOURIER_FOURIER` | `.FOURIER_FDTD` | `.SEM_SEM` | `.FOURIER_SEM`).
* `mainCreateConvergenceData.m`: Functionality for creating convergence plots similar to the experiments from the paper.
* `mainSingleMethods.m`: Implementation of the Fourier methods and SEM for running in the full domain (for investigation purposes).

* `mainAnalyzeErrors.m`: Functionality for assessing interface and overall errors similar to the experiments from the paper. 
* `mainAnalyzeLaplacianSideBySide.m`: Investigation of the level of interface errors occuring for calculating the Laplacian depending on the polynomial order in SEM.

* `result_plots/`: Various scripts for calculating efficiency (CPU) and other functionalities.