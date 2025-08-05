function [x, fval, iter] = projected_gradient(f, grad_f, proj, x0, max_iter, tol, verbose)
    
    % Implements the projected gradient descent with Armijo line search
    
    % Inputs:
    %   f      - Objective function handle
    %   grad_f - Gradient function handle
    %   proj   - Function handle for the projection operator.
    %   x0     - Initial guess
    %   max_iter - Maximum number of iterations
    %   tol    - Tolerance for convergence
    %   alpha  - Step size scaling factor (0 < alpha < 1)
    %   beta   - Step size reduction factor (0 < beta < 1)
    % Outputs:
    %   x      - Optimized variable
    %   fval   - Objective function value at x
    %   iter   - Number of iterations

    % Line search parameters
    sigma = 1e-4; % Armijo parameter
    beta = 0.25;   % Step size reduction factor

    % BB step-size parameters
    alphaBB = 1; 
    alphaBB_min = 1e-10;
    alphaBB_max = 1e4;
    alphaBB_M = 2;
    alphaBB2_hist = alphaBB * ones(max_iter, 1);
    tau = 0.1;

    x = x0;
    fval = f(x);
    g = grad_f(x);
    iter = 0;

    % Stagnation parameters
    max_iter_stag = 1000;
    count_stag = 0;
    best_norm = inf;
    best_x = [];
    best_fval = inf;

    while iter < max_iter

        d = proj(x - alphaBB * g) - x; % Compute projected gradient (descent direction)

        d_norm = norm(d);

        if d_norm < best_norm
            best_norm = d_norm;
            best_x = x;
            best_fval = fval;
            count_stag = 0;
        else
            count_stag = count_stag + 1;
        end

        if best_norm <= tol
            if verbose
                fprintf("\t\t Converged! \t minnorm(d): %.3f \n", best_norm);
            end
            x = best_x;
            fval = best_fval;
            break
        end

        if count_stag == max_iter_stag
            if verbose
                fprintf("\t\t Stagnation! \t minnorm(d): %.3f \n", best_norm);
            end
            x = best_x;
            fval = best_fval;
            break;
        end

        % Verbose logging
        if verbose && mod(iter, 100) == 0 && iter > 0
            fprintf('\t\t iter: %d \t norm(d): %.4f \t minnorm(d): %.4f\n', iter, d_norm, best_norm);
        end
        
        
        % Armijo line search
        alpha = 1;
        x_trial = x + alpha * d;
        f_trial = f(x_trial);
        while f_trial > fval + sigma * alpha * g' * d
            alpha = beta * alpha; % Reduce step size
            if alpha < 1e-20
                warning("Step size too small");
                break;
            end
            x_trial = x + alpha * d;
            f_trial = f(x_trial);
        end

        % Update variables
        x_new = x_trial;
        fval = f_trial;
        g_new = grad_f(x_new);

        % BB stepsize
        s = x_new - x;
        t = g_new - g;

        alphaBB1 = (s'*s)/(s'*t);
        alphaBB2 = (s'*t)/(t'*t);

        if iter > 1
            if s'*t < 0
                alphaBB = alphaBB_max;
            else
                alphaBB1 = max(alphaBB_min, min(alphaBB1, alphaBB_max));
                alphaBB2 = max(alphaBB_min, min(alphaBB2, alphaBB_max));
    
                if alphaBB2/alphaBB1 <= tau
                    j = max(1, iter - alphaBB_M);
                    alphaBB = min(alphaBB2_hist(j:iter));
                    tau = 0.9*tau;
                else
                    alphaBB = alphaBB1;
                    tau = 1.1*tau;
                end
            end
            alphaBB2_hist(iter) = alphaBB2;
        end

        x = x_new;
        g = g_new;
        iter = iter + 1;

    end
end