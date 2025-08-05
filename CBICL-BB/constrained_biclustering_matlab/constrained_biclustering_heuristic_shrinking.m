function [best_lb, best_Yu, best_Yv, time] = constrained_biclustering_heuristic_shrinking(Zshr, Tu, Tv, W, k, CL_u, CL_v, verbose, iter)

    [n, original_n] = size(Tu);
    [m, original_m] = size(Tv);
                
    Zuu = Zshr(1:n, 1:n);
    Zvv = Zshr(n+1:n+m, n+1:n+m); 
    
    [labels_u, ~] = kmeans(Zuu, k, 'Replicates', 500, 'Start', 'plus');
    Xu_bar = zeros(n, k);
    for i=1:n
        Xu_bar(i, labels_u(i)) = 1;
    end
    
    % force pairwise feasibility on Xu

    if iter == 0 % exact approach
        [Xu_f_shr, ~, timeU] = heuristic_one_way_pairwise_shrinking(Xu_bar, Tu, CL_u, verbose);
    else % heuristic approach
        [Xu_f_shr, ~, timeU] = greedy_assignment_CL(Xu_bar, Tu, CL_u, verbose, iter);
    end

    [labels_v, ~] = kmeans(Zvv, k, 'Replicates', 500, 'Start', 'plus');
    Xv_bar = zeros(m, k);
    for i=1:m
        Xv_bar(i, labels_v(i)) = 1;
    end

    % force pairwise feasibility on Xv

    if iter == 0 % exact approach
        [Xv_f_shr, ~, timeV] = heuristic_one_way_pairwise_shrinking(Xv_bar, Tv, CL_v, verbose);
    else % heuristic approach
        [Xv_f_shr, ~, timeV] = greedy_assignment_CL(Xv_bar, Tv, CL_v, verbose, iter);
    end

    time = timeU+timeV;

    if isempty(Xu_f_shr) || isempty(Xv_f_shr)

        best_Yu = [];
        best_Yv = [];
        best_lb = -inf;

    else

        Xu_f = Tu' * Xu_f_shr;
        Xv_f = Tv' * Xv_f_shr;

        %fprintf('\t Solving linear assignment problem...\n');

        % solve the linear assignment problem
    
        Yu_f = Xu_f*diag(1./sqrt(sum(Xu_f, 1)));
        Yv_f = Xv_f*diag(1./sqrt(sum(Xv_f, 1)));
        
        % Index helper function
        assidx = @(i, j) i+(j-1)*k;
        
        % Build model
        model.A = sparse(2*k, k*k);
        model.lb = zeros(k*k, 1);
        obj = zeros(k, k);
    
        for i=1:k
            for j=1:k
                model.A(i, assidx(i, j)) = 1;
                obj(i, j) = Yu_f(:, i)'*W*Yv_f(:, j);
            end
        end
            
            
        for j=1:k
            model.A(j+k, 1+k*(j-1):k*j) = 1;
        end
            
        model.rhs = ones(2*k, 1);
        model.sense = repmat('=', 2*k, 1);
                
        model.modelsense = 'max';
        model.obj = obj(:);
        params.outputflag = verbose;
        result = gurobi(model, params);
        %disp(result.objval);
        delta = reshape(result.x, [k, k]);
        [r_idx, c_idx] = find(delta);
            
        % recover best normalized partition matrices
        best_Yu = zeros(original_n, k);
        best_Yv = zeros(original_m, k);
        for i=1:k
            best_Yu(:, i) = Yu_f(:, r_idx(i));
            best_Yv(:, i) = Yv_f(:, c_idx(i));   
        end
        
        best_lb = trace(best_Yu'*W*best_Yv);

    end
     
end