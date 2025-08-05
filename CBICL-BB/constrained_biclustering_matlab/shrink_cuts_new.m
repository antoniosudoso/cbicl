function new_Bcell = shrink_cuts_new(ML_graph_U, ML_graph_V, init_B_cell, n, m, original_n)
    
    n_ineq = size(init_B_cell, 2);
    
    new_Bcell = cell(1, n_ineq);
    [bins_v_U , ~] = conncomp(ML_graph_U);
    [bins_v_V , ~] = conncomp(ML_graph_V);

    counter = 1;
    for c=1:n_ineq
        
            [id_i, id_j, v] = find(init_B_cell{c});
                
            d = size(id_i, 1);

            for t=1:d
                if id_i(t) <= original_n
                    id_i(t) = bins_v_U(id_i(t));
                else
                    id_i(t) = bins_v_V(id_i(t) - original_n) + n;
                end
                if id_j(t) <= original_n
                    id_j(t) = bins_v_U(id_j(t));
                else
                    id_j(t) = bins_v_V(id_j(t) - original_n) + n;
                end
            end
            new_Bcell{counter} = sparse(id_i, id_j, v, n + m, n + m);
            counter = counter + 1;

    end
    
    n_ineq = counter - 1;
    new_Bcell = new_Bcell(1:n_ineq);
    
end