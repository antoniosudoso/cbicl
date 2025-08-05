clear;

addpath(genpath('/home/antonio/Dropbox/biclustering/SDPNAL+/'));
addpath(genpath('/home/antonio/gurobi1100'));

params = struct();
params.n_threads = 20;
params.bb_tol = 1e-3;
params.sdp_tol = 1e-4;
params.sdp_verbose = 0;
params.gurobi_verbose = 0;
params.heuristic_iter = 0;
params.cp_maxiter = 10;
params.cp_tol = 1e-4;
params.cp_epsineq = 1e-4;
params.cp_maxineq = 100000;
params.cp_percineq = 0.10;
params.cp_activeineq = 1e-4;

A = readmatrix("/home/antonio/Dropbox/biclustering/DeSouto/Bhattacharjee-2001-v2.txt");
[ML_U, ML_V, CL_U, CL_V] = parse_constraints("/home/antonio/Dropbox/biclustering/DeSouto/constraints_sdp/2/seed_0/Bhattacharjee-2001-v2_4_MLu_CLu_MLv_CLv.txt");
k = 4;

result = call_solve_constrained_biclustering_root(A, k, ML_U, CL_U, ML_V, CL_V, params);
