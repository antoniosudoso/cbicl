function is_feasible = check_feasibility(Xu, Xv, ML_u, CL_u, ML_v, CL_v)

    % Precompute inner products
    Xuu = Xu * Xu';
    Xvv = Xv * Xv';

    % Helper to compute percentage satisfied constraints
    function percent = check_constraint(mat, pairs, condition)
        if isempty(pairs)
            percent = 100;
            return;
        end
        vals = mat(sub2ind(size(mat), pairs(:,1), pairs(:,2)));
        if strcmp(condition, '>0')
            satisfied = vals > 0;
        elseif strcmp(condition, '=0')
            satisfied = vals == 0;
        else
            error('Unsupported condition');
        end
        percent = 100 * sum(satisfied) / size(pairs, 1);
    end

    % Check and report
    ml_u_pct = check_constraint(Xuu, ML_u, '>0');
    cl_u_pct = check_constraint(Xuu, CL_u, '=0');
    ml_v_pct = check_constraint(Xvv, ML_v, '>0');
    cl_v_pct = check_constraint(Xvv, CL_v, '=0');

    fprintf('\t ML_u satisfied = %.2f (%%)\n', ml_u_pct);
    fprintf('\t CL_u satisfied = %.2f (%%)\n', cl_u_pct);
    fprintf('\t ML_v satisfied = %.2f (%%)\n', ml_v_pct);
    fprintf('\t CL_v satisfied = %.2f (%%)\n', cl_v_pct);

    % Return true if all constraints are fully satisfied
    is_feasible = (ml_u_pct == 100) && (cl_u_pct == 100) && ...
                  (ml_v_pct == 100) && (cl_v_pct == 100);
end
