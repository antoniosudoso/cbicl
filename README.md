## Exact and Heuristic Algorithms for Constrained Biclustering

This repository contains the implementation of **CBICL-BB** and **CBICL-LR**, described in the paper ["Exact and Heuristic Algorithms for Constrained Biclustering"](https://arxiv.org/abs/2508.05493). **CBICL-BB** is an exact algorithm, based on the branch-and-cut technique, for constrained biclustering through the $k$-densest-disjoint biclique criterion with must-link and cannot-link constraints. **CBICL-LR** is a heuristic algorithm based on a low-rank factorization of the SDP relaxation.  Datasets and constraint sets used in the experiments are available in the `Data` folder.

## CBICL-BB

### Installation
**CBICL-BB** is implemented in C++ with some routines written in MATLAB. **CBICL-BB** calls [SDPNAL+](https://blog.nus.edu.sg/mattohkc/softwares/sdpnalplus/) for solving SDP relaxations and [Gurobi](https://www.gurobi.com/) as integer linear programming solver for the rounding heuristic.


Ubuntu and Debian instructions:

1) Install MATLAB (>= 2021b)

2) Install Gurobi (>= 10.0.2)

3) Install CMake, OpenBLAS, LAPACK and Armadillo:
 ```
sudo apt-get update
sudo apt-get install cmake libopenblas-dev liblapack-dev libarmadillo-dev
```
4) Open the makefile `CBICL-BB/constrained_biclustering_cpp/Makefile` 
	- Set the variable `matlab_path` with your MATLAB installation folder.

5) Compile the code:

```
cd CBICL-BB/constrained_biclustering_cpp/
make
```

4) Download SDPNAL+, move the folder `CBICL-BB/constrained_biclustering_matlab` (containing the MATLAB source code of **CBICL-BB**) into the SDPNAL+ main directory and set the parameter `SDP_SOLVER_FOLDER` in the configuration file accordingly. This folder and its subfolders will be automatically added to the MATLAB search path when the algorithm starts.

This code has been tested under Ubuntu 22.04 LTS with MATLAB R2022b, Gurobi 11.0.0

### Configuration
Various parameters used in **CBICL-BB** can be modified in the configuration file `CBICL-BB/constrained_biclustering_cpp/config.txt`:

- `BRANCH_AND_BOUND_TOL` - optimality tolerance
- `BRANCH_AND_BOUND_PARALLEL` -  thread pool size: single thread (1), multi-thread (> 1)
- `BRANCH_AND_BOUND_MAX_NODES` - maximum number of nodes
- `BRANCH_AND_BOUND_VISITING_STRATEGY` - best first (0),  depth first (1), breadth first (2)
- `MATLAB_SESSION_THREADS_ROOT` - number of threads for the MATLAB session at the root noee
- `MATLAB_SESSION_THREADS_CHILD` - number of threads for the MATLAB session for child nodes
- `SDP_SOLVER_FOLDER` - full path of SDPNAL+ folder
- `SDP_SOLVER_TOL` - accuracy of SDPNAL+ in the relative KKT residual
- `SDP_SOLVER_VERBOSE` - do not display log (0), display log (1)
- `CP_MAX_ITER` - maximum number of cutting-plane iterations
- `CP_TOL` - tolerance between two consecutive cutting-plane iterations
- `CP_MAX_INEQ` - maximum number of valid inequalities to separate
- `CP_PERC_INEQ` - fraction of the most violated inequalities to add
- `CP_EPS_INEQ` - tolerance for checking the violation of the inequalities
- `CP_EPS_ACTIVE` - tolerance for detecting active inequalities
- `GUROBI_FOLDER` - Gurobi solver path
- `GUROBI_VERBOSE` - do not display log (0), display log (1)

### Usage
```
cd biclustering_cpp/
./bb <W_PATH> <K> <CONSTRAINTS_PATH> <LOG_PATH> <RESULT_PATH>
```
- `W_PATH` - data matrix
- `K` - number of biclusters
- `CONSTRAINTS_PATH` - pairwise constraints
- `LOG_PATH` - log file
- `RESULT_PATH` - optimal assignment matrices

File `W_PATH` contains the weights `w_ij` and the must include a header line with the number of rows `n` and columns `m`:

```
n m
w_11 w_12 ... w_1m
w_21 w_22 ... w_2m
...
...
w_n1 w_n2 ... w_nm
```

File `CONSTRAINTS_PATH` includes indices `(i, j)` of the entities involved in must-link (ML_U/ML_V) and/or cannot-link (CL_U/CL_V) constraints:

```
CL_U i1 j1
CL_V i2 j2
...
...
ML_U i3 j3
ML_V i4 j4
```

Results are saved in `RESULT_PATH` producing two matrix, the `n x k` and `m x k` row and column assignments, respectively.

### Log

The log file reports the progress of the algorithm:

- `N` - number of rows at the current node
- `M` - number of columns at the current node
- `ID_PAR` - id of the parent node
- `ID` - id of the current node
- `UB_PAR` - upper bound of the parent node
- `UB` - upper bound of the current node
- `TIME (s)` - running time in seconds of the current node
- `CP_ITER` - number of cutting-plane iterations
- `CP_FLAG` - termination flag of the cutting-plane procedure
    - `-2` - maximum number of iterations
    - `-1` - SDP not solved or partially solved successfully
    -  ` 0` - no violated inequalities
    -  ` 1` - node must be pruned
    -  ` 2` - upper bound greater than the previous one
    -  ` 3` - upper bound decrease is not sufficiently large
- `CP_INEQ` - number of inequalities added in the last cutting-plane iteration
- `LB` - current lower bound
- `HTIME (s)` - running time in seconds of the rounding heuristic
- `BEST_LB` - global lower bound
- `SET` - vertex set selection for branching
    -  `U` - branch on the vertices in U
    -  `V` - branch on the vertices in V
    -  `-1` - branching is not needed
- `I J` - indices of branching decision
- `NODE_GAP` - gap at the current node
- `GAP` - overall gap 
- `OPEN` - number of open nodes


## CBICL-LR

### Usage
**CBICL-BB** is fully implemented in MATLAB and can be run by using:

```
[YU_out, YV_out, result] = bm_constrained_biclustering(A, k, r, ML_U, CL_U, ML_V, CL_V, params, YU_init, YV_init)
```
where:

- `A` - `n x m` data matrix
- `k` - number of biclusters
- `r` - target rank (automatically chosen if empty [])
- `ML_U` - `|ML_U| x 2` matrix of must-link constraint for the vertex set U
- `CL_U` - `|CL_U| x 2` matrix of cannot-link constraint for the vertex set U
- `ML_V` - `|ML_V| x 2` matrix of must-link constraint for the vertex set V
- `CL_V` - `|CL_V| x 2` matrix of cannot-link constraint for the vertex set V
- `params` - struct of parameters for the augmented Lagrangian method
	- `params.tol` - tolerance of stopping criteria
	- `params.primal_tol` - tolerance for the subproblem
	- `params.maxiter` - maximum of iterations
	- `params.tau` - parameter of sufficient decrease for updating beta
	- `params.gamma` - parameter for updating beta
	- `params.maxtime` - time limit in seconds
- `YU_init` - initial `n x r` solution for the vertex set U (randomly chosen if empty [])
- `YV_init` - initial `m x r` solution for the vertex set V (randomly chose if empty [])

The output is given by:

- `YU_out` - `n x r` solution for the vertex set U
- `YV_out` - `m x r` solution for the vertex set V
- `result` - struct of computed metrics
	- `result.r` - chosen rank
	- `result.inner_iter` - number of projected gradient iterations
    - `result.outer_iter` - number of augmented Lagrangian iterations
    - `result.infeas_list` - residuals history
    - `result.beta_list` - beta history
    - `result.time` - computational time in seconds
