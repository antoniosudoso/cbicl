function [X_pairwise, obj, time] = heuristic_one_way_pairwise(Xbar, ML, CL, verbose)

    tStart = tic;

    % Xbar: unconstrained row or column assignment matrix n x k or m x k
    % ML: must-link constraints matrix n_ml x 2
    % CL: cannot-link constraints matrix n_cl x 2
    % verbose: display or not gurobi log
    
    [n, k] = size(Xbar);
    
    n_ml = size(ML, 1);
    n_cl = size(CL, 1);
    
    % Index helper function
    assidx = @(i, j) i+(j-1)*n;
    
    % Build model
    model.vtype = repmat('B', n*k, 1);
    model.A = sparse(n+k*n_ml+k*n_cl, n*k);
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
    for c=1:n_ml
        i = ML(c, 1);
        j = ML(c, 2);
        for h=1:k
            model.A(count, assidx(i, h)) = 1;
            model.A(count, assidx(j, h)) = -1;
            count = count + 1;
        end
    end
    for c=1:n_cl
        i = CL(c, 1);
        j = CL(c, 2);
        for h=1:k
            model.A(count, assidx(i, h)) = 1;
            model.A(count, assidx(j, h)) = 1;
            count = count + 1;
        end
    end
    % spy(model.A)
    model.rhs = [ones(n, 1); ones(k, 1); zeros(k*n_ml, 1); ones(k*n_cl, 1)]; % n + k + k*n_ml + k*n_cl constraints
    model.sense = [repmat('=', n, 1); repmat('>', k, 1); repmat('=', k*n_ml, 1); repmat('<', k*n_cl, 1)];
        
    %fprintf('Finding nearest assignment matrix...\n');
    model.modelsense = 'max';
    model.obj = Xbar(:);
    params.outputflag = verbose;
    % params.threads = n_threads;
    result = gurobi(model, params);
    % disp(result.objval);
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