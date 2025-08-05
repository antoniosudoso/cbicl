function result = call_solve_constrained_biclustering_root(W, k, ML_u, CL_u, ML_v, CL_v, params)

    %disp(params)

    fprintf("\n\t Number of user ML_u constraints: %d \n", size(ML_u, 1));
    fprintf("\t Number of user CL_u constraints: %d \n", size(CL_u, 1));
    fprintf("\t Number of user ML_v constraints: %d \n", size(ML_v, 1));
    fprintf("\t Number of user CL_v constraints: %d \n", size(CL_v, 1));

    [n, m] = size(W);

    result = solve_constrained_biclustering_shrinking(W, k, ML_u, CL_u, ML_v, CL_v, cell(0), -inf, [], [], params);
    result.htime = sum(result.htime_list);


    if result.cp_flag == 1
        result.branching_type = -1;
        result.i_idx = -1;
        result.j_idx = -1;
    else
        [branching_type, i_idx, j_idx, ~] = get_branching_pair(result.best_Z, n, m);
        result.branching_type = branching_type;
        result.i_idx = i_idx;
        result.j_idx = j_idx;
    end
    
    %disp(result)
    
end