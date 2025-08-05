function [lb, Yu, Yv, tot_time, avg_time, avg_iter, avg_infeas, avg_time_ILP] = bm_multistart(n_times, A, k, ML_U, CL_U, ML_V, CL_V, params)

    tot_time = 0;
    lb = -inf;
    iter_list = [];
    infeas_list = [];
    time_ILP_list = [];
    Yu = [];
    Yv = [];

    for i=1:n_times
    
        [YU_out, YV_out, result] = bm_constrained_biclustering(A, k, [], ML_U, CL_U, ML_V, CL_V, params, [], []);
        tot_time = tot_time + result.time;
        iter_list = [iter_list; result.outer_iter];
        infeas_list = [infeas_list; result.infeas];
    
        Z = [YU_out*YU_out', YU_out*YV_out'; YV_out*YU_out', YV_out*YV_out'];
        [best_lb, YU, YV, timeILP] = constrained_biclustering_heuristic_one_way(Z, A, k, ML_U, CL_U, ML_V, CL_V, 0, 1);
        if best_lb > lb
            lb = best_lb;
            Yu = YU;
            Yv = YV;
        end


        time_ILP_list = [time_ILP_list; timeILP];

    end

    avg_time = tot_time / n_times;
    avg_iter = mean(iter_list);
    avg_infeas = mean(infeas_list);
    avg_time_ILP = mean(time_ILP_list);

end
