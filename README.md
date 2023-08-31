# Numerical PDE solvers
Author: Nikolas Borrel-Jensen (2023).

The repository consists of three projects:

* `SEMSolvers/`: Spectal-element solvers in 1D and 2D for solving the acoustic wave equation and the Helmholtz equations with frequency-independent and dependent boundary conditions.
* `FourierSEMSolver/`: The Fourier-SEM domain decomposition method for coupling the Fourier method with the spectral elemement method through a finite-difference time-domain scheme.
* `WBBEMSolver/`: Wave-based boundary element solver (also known as the 'WBM') for solving the Helmholtz equation in 2D.
* `shared/`: Contains shared functionality and should be added to the MATLAB search path `addpath('shared')`.
* `scripts/`: Scripts for running on the DTU HPC cluster and other `json` setup files.

A README describing each project can be found in the root of the folders.