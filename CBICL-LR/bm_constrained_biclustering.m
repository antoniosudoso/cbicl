function [YU_out, YV_out, result] = bm_constrained_biclustering(A, k, r, init_ML_u, init_CL_u, init_ML_v, init_CL_v, params, init_YU, init_YV)

    %rng(1993)

    [original_n, original_m] = size(A); % Number of rows and columns
    At = A';

    % transformation for shrinking U
    [TU, ML_graph_U] = build_T(original_n, init_ML_u);
    n = size(TU, 1); % new size
    CL_U = update_CL(ML_graph_U, init_CL_u);
    n_CL_U = size(CL_U, 1);

    TUt = TU';
    e_n = ones(n, 1);
    eU = TU*ones(original_n, 1);
    TUTUt = TU*TUt;
    
    % transformation for shrinking V
    [TV, ML_graph_V] = build_T(original_m, init_ML_v);
    m = size(TV, 1); % new size
    CL_V = update_CL(ML_graph_V, init_CL_v);
    n_CL_V = size(CL_V, 1);

    TVt = TV';
    e_m = ones(m, 1);
    eV = TV*ones(original_m, 1);
    TVTVt = TV*TVt;

    TVAtTUt = TV*At*TUt;
    TUATVt = TU*A*TVt;
    TUATVt_vec = TUATVt(:);


    if isempty(r)
        n_constr = n+1+n_CL_U+m+1+n_CL_V;
        r = ceil((-1 + sqrt(1 + 8*n_constr)) / 2);
    end

    result.r = r;

    vec = @(X) reshape(X, [], 1);
    matU = @(x) reshape(x, n, r);
    matV = @(x) reshape(x, m, r);
    proj = @(x) max(0, x);    

    CL_U_op = @(X) X(sub2ind([n, n], CL_U(:, 1), CL_U(:, 2)));
    CL_U_E = @(w) sparse(CL_U(:, 1), CL_U(:, 2), w, n, n);
    CL_V_op = @(X) X(sub2ind([m, m], CL_V(:, 1), CL_V(:, 2)));
    CL_V_E = @(w) sparse(CL_V(:, 1), CL_V(:, 2), w, m, m);

    beta = 10; % Coeffient of the augmented term

    if isempty(init_YU)
        YU_0 = abs(randn(n, r)); 
        YU = proj(YU_0);
    else
        if size(init_YU, 2) == k
            % Fill with zeros
            YU = zeros(n, r);
            YU(:, 1:k) = inv(TUTUt)*TU*init_YU;
        else
            YU = init_YU;
        end
    end

    if isempty(init_YV)
        YV_0 = abs(randn(m, r)); 
        YV = proj(YV_0);
    else
        if size(init_YV, 2) == k
            % Fill with zeros
            YV = zeros(m, r);
            YV(:, 1:k) = inv(TVTVt)*TV*init_YV;
        else
            YV = init_YV;
        end
    end

    yU = vec(YU);
    yV = vec(YV);
    lU = zeros(n, 1);
    lV = zeros(m, 1);
    mU = 0;
    mV = 0;
    wU = zeros(n_CL_U, 1);
    wV = zeros(n_CL_V, 1);

    lU_min = -1e20*ones(n, 1);
    lU_max =  1e20*ones(n, 1);
    lV_min = -1e20*ones(m, 1);
    lV_max =  1e20*ones(m, 1);
    mU_min = -1e20;
    mU_max =  1e20;
    mV_min = -1e20;
    mV_max =  1e20;
    wU_min = -1e20*ones(n_CL_U, 1);
    wU_max =  1e20*ones(n_CL_U, 1);
    wV_min = -1e20*ones(n_CL_V, 1);
    wV_max =  1e20*ones(n_CL_V, 1);

    safeguard = @(x, lb, ub) min(max(lb, x), ub);

    log_proj = 0;

    infeasU_sum = matU(yU)*(matU(yU)'*eU)-e_n;
    infeasV_sum = matV(yV)*(matV(yV)'*eV)-e_m;
    infeasU_norm = trace(TUTUt*matU(yU)*matU(yU)')-k;
    infeasV_norm = trace(TVTVt*matV(yV)*matV(yV)')-k;
    infeasU_CL = CL_U_op(matU(yU)*matU(yU)');
    infeasV_CL = CL_V_op(matV(yV)*matV(yV)');
    infeas = max(norm([infeasU_sum; infeasU_norm; infeasU_CL])/sqrt(n+1+n_CL_U), norm([infeasV_sum; infeasV_norm; infeasV_CL])/sqrt(m+1+n_CL_U));
    
    fprintf("\n");
    
    result.inner_iter = 0;
    result.outer_iter = 0;
    result.infeas_list = [];
    result.beta_list = [];

    tStart = tic;

    %keyboard


    for i=1:params.maxiter

        result.inner_iter = result.inner_iter + 1;

        %fprintf("\n\t------------------------------- Updating yU ----------------------------------\n")
        ffU_vec = @(y_U) augL_constrained(y_U, yV, lU, lV, mU, mV, wU, wV, beta, TUATVt_vec, k, matU, matV, TUTUt, TVTVt, eU, eV, e_n, e_m, CL_U_op, CL_V_op);
        fgradU_vec = @(y_U) vec(gradU_constrained(y_U, yV, lU, mU, wU, beta, TUATVt, k, matU, matV, TUTUt, eU, e_n, CL_U_op, CL_U_E));
        [yU_new, ~, ~] = projected_gradient(ffU_vec, fgradU_vec, proj, yU, params.maxiter, params.primal_tol, log_proj);
        %fprintf("\n\t------------------------------- Updating yV ----------------------------------\n")
        ffV_vec = @(y_V) augL_constrained(yU_new, y_V, lU, lV, mU, mV, wU, wV, beta, TUATVt_vec, k, matU, matV, TUTUt, TVTVt, eU, eV, e_n, e_m, CL_U_op, CL_V_op);
        fgradV_vec = @(y_V) vec(gradV_constrained(yU_new, y_V, lV, mV, wV, beta, TVAtTUt, k, matU, matV, TVTVt, eV, e_m, CL_V_op, CL_V_E));
        [yV_new, ~, ~] = projected_gradient(ffV_vec, fgradV_vec, proj, yV, params.maxiter, params.primal_tol, log_proj);

        rdiffU = norm(yU_new - yU)/sqrt(n*r);
        rdiffV = norm(yV_new - yV)/sqrt(m*r);
        rdiff = max(rdiffU, rdiffV);

        %fprintf("\t Primal update with beta %d \t rdiffU %.4f \t rdiffV %.4f \t rdiff %.4f \n", beta, rdiffU, rdiffV, rdiff);


        if rdiff <= params.tol

            result.outer_iter = result.outer_iter + 1;
            
            infeasU_sum_new = matU(yU_new)*(matU(yU_new)'*eU)-e_n;
            infeasV_sum_new = matV(yV_new)*(matV(yV_new)'*eV)-e_m;
            infeasU_norm_new = trace(TUTUt*matU(yU_new)*matU(yU_new)')-k;
            infeasV_norm_new = trace(TVTVt*matV(yV_new)*matV(yV_new)')-k;
            infeasU_CL_new = CL_U_op(matU(yU_new)*matU(yU_new)');
            infeasV_CL_new = CL_V_op(matV(yV_new)*matV(yV_new)');
            infeas_new = max(norm([infeasU_sum_new; infeasU_norm_new; infeasU_CL_new])/sqrt(n+1+n_CL_U), norm([infeasV_sum_new; infeasV_norm_new; infeasV_CL_new])/sqrt(m+1+n_CL_V));

            result.infeas_list = [result.infeas_list; infeas_new];
            result.beta_list = [result.beta_list; beta];

            fprintf("\n\t Dual update with beta %.4e \t infeasU %.4f \t infeasV %.4f \t infeas %.4f ", beta, ...
                norm([infeasU_sum_new; infeasU_norm_new; infeasU_CL_new])/sqrt(n+1+n_CL_U), ...
                norm([infeasV_sum_new; infeasV_norm_new; infeasV_CL_new])/sqrt(m+1+n_CL_V), infeas_new);

            % Stopping criteria    
            if infeas_new <= params.tol
               fprintf("\n\t Converged! \n\n");
               result.infeas = infeas_new;
               break
            end

            % Dual update
            lU = safeguard(lU+beta*infeasU_sum_new, lU_min, lU_max);
            lV = safeguard(lV+beta*infeasV_sum_new, lV_min, lV_max);
            mU = safeguard(mU+beta*infeasU_norm_new, mU_min, mU_max);
            mV = safeguard(mV+beta*infeasV_norm_new, mV_min, mV_max);
            wU = safeguard(wU+beta*infeasU_CL_new, wU_min, wU_max);
            wV = safeguard(wV+beta*infeasV_CL_new, wV_min, wV_max);

            %fprintf("\t infeas_new %.4f \t infeas_old %.4f \n", infeas_new, infeas);
            % Update penalty parameter 
            if infeas_new > params.tau * infeas
                beta = min(params.gamma * beta, 1e14);
            end

            infeas = infeas_new;
            result.infeas = infeas;


        end

        yU = yU_new;
        yV = yV_new;


        if toc(tStart) > params.maxtime
            fprintf("\n\t Time limit reached!\n\n");
            break
        end


    end

    result.time = toc(tStart); 

    % Output solution
    YU_out = TUt*matU(yU);
    YV_out = TVt*matV(yV);


end