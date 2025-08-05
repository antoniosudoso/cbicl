function [X_pairwise, obj, time] = heuristic_one_way_pairwise_shrinking(Xbar_shr, T, CL, verbose)

    tStart = tic;

    % Xbar_shr: unconstrained row or column assignment matrix n x k or m x k
    % CL: cannot-link constraints matrix n_cl x 2
    % T: transformation matrix
    % verbose: display or not gurobi log
    
    [n, k] = size(Xbar_shr);
    
    n_cl = size(CL, 1);
    
    % Index helper function
    assidx = @(i, j) i+(j-1)*n;
    
    % Build model
    model.vtype = repmat('B', n*k, 1);
    model.A = sparse(n+k+k*n_cl, n*k);
    count = 1;
    for i=1:n
        for j=1:k
            model.A(i, assidx(i, j)) = 1;
        end
        count = count + 1;
    end
    for j=1:k
        model.A(count, 1+n*(j-1):n*j) = 1;
        count = count + 1;
    end

    for c=1:n_cl
        i = find(T(:, CL(c, 1)));
        j = find(T(:, CL(c, 2)));
        for h=1:k
            model.A(count, assidx(i, h)) = 1;
            model.A(count, assidx(j, h)) = 1;
            count = count + 1;
        end
    end
    % spy(model.A)
    model.rhs = [ones(n, 1); ones(k, 1); ones(k*n_cl, 1)]; % n + k + k*n_cl constraints
    model.sense = [repmat('=', n, 1); repmat('>', k, 1); repmat('<', k*n_cl, 1)];
        
    %fprintf('Finding nearest assignment matrix...\n');
    model.modelsense = 'max';
    model.obj = Xbar_shr(:);
    params.outputflag = verbose;
    % params.threads = n_threads;
    result = gurobi(model, params);
    if result.status ~= "OPTIMAL"
        X_pairwise = [];
        obj = -inf;
    else
        X_pairwise = reshape(result.x, [n, k]);
        obj = result.objval;
    end

    tEnd = toc(tStart);
    time = tEnd;
     
end