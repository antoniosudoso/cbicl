function result = call_solve_constrained_biclustering_child(W, k, ML_u, CL_u, ML_v, CL_v, init_B_cell, global_lb, global_Xu, global_Xv, params)
    
    %disp("ML_u")
    %disp(ML_u)
    %disp("CL_u")
    %disp(CL_u)
    %disp("ML_v")
    %disp(ML_v)
    %disp("CL_v")
    %disp(CL_v)
    
    %fprintf("\n\t Number of child ML_u constraints: %d \n", size(ML_u, 1));
    %fprintf("\t Number of child CL_u constraints: %d \n", size(CL_u, 1));
    %fprintf("\t Number of child ML_v constraints: %d \n", size(ML_v, 1));
    %fprintf("\t Number of child CL_v constraints: %d \n", size(CL_v, 1));

    [n, m] = size(W);

%     is_feasible_u = is_feasible_one_way_pairwise(n, k, ML_u, CL_u, params.gurobi_verbose);
%     is_feasible_v = is_feasible_one_way_pairwise(m, k, ML_v, CL_v, params.gurobi_verbose);
% 
%     if is_feasible_u && is_feasible_v
% 
%         result = solve_constrained_biclustering_shrinking(W, k, ML_u, CL_u, ML_v, CL_v, init_B_cell, global_lb, global_Xu, global_Xv, params);
%         result.best_lb = max(result.lb_list);
% 
%     else
%         
%         result.best_lb = -inf;
%         result.best_ub = inf;
%         result.best_Xu = [];
%         result.best_Xv = [];
%         result.cp_iter = 0;
%         result.cp_flag = 1;
%         result.best_B_cell = cell(0);
%         result.n = n;
%         result.m = m;
% 
%     end

    result = solve_constrained_biclustering_shrinking(W, k, ML_u, CL_u, ML_v, CL_v, init_B_cell, global_lb, global_Xu, global_Xv, params);
    result.best_lb = max(result.lb_list);
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