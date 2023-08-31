# Spectral-Element solvers
Code written by Nikolas Borrel-Jensen adjusted from existing code: the spectral element helper routines (inside `shared/`) are mainly written by Allan P. Engsig-Karup and the original code for handling frequency-dependent boundary conditions for the SEM is written by Finnur Pind.

## Wave equation solvers
The main codes for calculating solutions to the wave equation in 1D and 2D are located here:

* `SEM_WEQ1D.m`: 1D solver
* `SEM_WEQ2D.m`: 2D solver

The solvers can be used for generating data and has been used for the 1D and 2D examples for the following papers:

* [Sound propagation in realistic interactive {3D} scenes with parameterized sources using deep neural operators](https://arxiv.org/abs/2308.05141) (N. Borrel-Jensen, S. Goswami, A. P. Engsig-Karup, G.E. Karniadakis, C.-H. Jeong, Arxiv 2023)
* [A sensitivity analysis on the effect of hyperparameters in deep neural operators applied to sound propagation](https://www.fa2023.org) (N. Borrel-Jensen, A. P. Engsig-Karup, C.-H. Jeong, 10th Convention of the European Acoustic Association of Forum Acusticum 2023)

### Generating 2D data for DeepONet model training 
The most important settings for generating data are explained below:

**Environment**
* `hpc_env`: set to `true` to run on DTU HPC cluster
* `run_parallel`: set to `true` to run in parallel. The line `parfor i=1:size(p_ics,1)` also needs to be set accordingly.
* `sim_id`: the id of the simulation output.

**Initial and boundary conditions**
* `boundary_type`: `freq_indep` | `neumann` | `dirichlet` | `freq_indep` | `freq_dep`
* `source_type`: can be a Gaussian pulse (used for all experiments) or a Gaussian Random Field.

**Geometry**
* `geometry`: three geometries are supported; L-shape, rectangular domain, 'furnished' rectangular domain

**Geometry sizes (branch net)**
* `xminmax = [0.0,3.0]`: trunk net input dimension (non-uniform grid) for x dimension. For non-rectangular domains, the size is determining the outer dimension and is 
adjusted in code.
* `yminmax = [0.0,3.0]`: same as above but for y dimension.

**Geometry sizes (trunk net)**
* `xminmax_u = [0.0,3.0]`: branch net input dimension (uniform grid) for x dimension. Could differ from trunk net input when used for transfer learning where the target model should have same branch input size as the source model.
* `yminmax_u = [0.0,3.0]`: same as above but for y dimension.

**Resolutions**
* `ppw_u`: points per wavelength for branch net input function `u` (default 2).
* `ppw`: points per wavelength for the spatial resolution as input to the trunk net (default 4-6).
* `dt_fixed`: the temporal resolution can be hardcoded and is used to ensure same temporal resolutions for training and validation/test data. I.e., in case the spatial resolution of the training data is finer than the validation/test data,  `dt_fixed` should be set equal to the training data resolution and the Courant condition will still be satisfied. If -1 is set, the property is ignored.

**Sources**
* `fixed_src_pos`: boolean determining if the source positions should be hardcoded or determined from source densities.
* `src_pad`: The Gaussian pulse has a width and therefore a threshold boundary padding is required for the Gaussian shape to be well-defined.
* `max_num_srcs`: when the area where the source position is allowed to move freely is narrow, the meshing algorithm might distribute the points mainly on the boundaries. Hence, it can be necessary to sample more densely and afterwards randomly extract source locations corresponding to the required source densitiy. If -1 is set, the property is ignored.

**Final notes**
* The data is generated with the resolution defined in the setup and required by the Courant condition needed for exact and stable simulations. However, the consumer of the data (e.g., for scientific machine learning) might required coarser resolution and therefore Python scripts have been made to prune the data and convert from 32 bit float to 16 bit float. The data can also be pruned directly from this code by setting `ppw_x_out` and `ppw_t_out`.

## Helmholtz solvers

The main codes for calculating solutions to the Helmholtz equation in 1D and 2D are located here:

* `SEM_Helmholtz1D.m`: 1D solver
* `SEM_Helmholtz2D.m`: 2D solver