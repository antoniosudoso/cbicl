function [best_lb, best_Yu, best_Yv, time_ILP] = constrained_biclustering_heuristic_one_way(Z, W, k, ML_u, CL_u, ML_v, CL_v, verbose, use_Zuv)

    [n, m] = size(W);
                
    Zuu = Z(1:n, 1:n);
    Zvv = Z(n+1:n+m, n+1:n+m); 
    Zuv = Z(1:n, n+1:n+m);
    
    if use_Zuv
        [labels_u, ~] = kmeans(Zuv, k, 'Replicates', 100, 'Start', 'plus');
    else
        [labels_u, ~] = kmeans(Zuu, k, 'Replicates', 100, 'Start', 'plus');
    end
    Xu_bar = zeros(n, k);
    for i=1:n
        Xu_bar(i, labels_u(i)) = 1;
    end

    tStart = tic;

    % force pairwise feasibility on Xu
    [Xu_f] = heuristic_one_way_pairwise(Xu_bar, ML_u, CL_u, 0);

    tEnd1 = toc(tStart);


    if use_Zuv
        [labels_v, ~] = kmeans(Zuv', k, 'Replicates', 100, 'Start', 'plus');
    else
        [labels_v, ~] = kmeans(Zvv, k, 'Replicates', 100, 'Start', 'plus');
    end
    Xv_bar = zeros(m, k);
    for i=1:m
        Xv_bar(i, labels_v(i)) = 1;
    end


    tStart = tic;

    % force pairwise feasibility on Xv
    [Xv_f] = heuristic_one_way_pairwise(Xv_bar, ML_v, CL_v, 0);

    tEnd2 = toc(tStart);

    time_ILP = tEnd1 + tEnd2;


    if isempty(Xu_f) || isempty(Xv_f)

        best_Yu = [];
        best_Yv = [];
        best_lb = -inf;

    else

        fprintf('\t Solving linear assignment problem...\n');

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
        best_Yu = zeros(n, k);
        best_Yv = zeros(m, k);
        for i=1:k
            best_Yu(:, i) = Yu_f(:, r_idx(i));
            best_Yv(:, i) = Yv_f(:, c_idx(i));   
        end
        
        best_lb = trace(best_Yu'*W*best_Yv);

    end
     
end