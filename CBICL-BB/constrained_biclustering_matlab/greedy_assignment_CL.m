function [X, f, time] = greedy_assignment_CL(Xbar_shr, T, CL_or, verbose, numRuns)

    % Adjust CL_or according to the connected components 
    n_cl = size(CL_or, 1);
    CL = zeros(n_cl, 2);
    for c=1:n_cl
        i = find(T(:, CL_or(c, 1)));
        j = find(T(:, CL_or(c, 2)));
        CL(c, 1) = i;
        CL(c, 2) = j;
    end

    tStart = tic;

    [n, k] = size(Xbar_shr);
    X_best = zeros(n, k);
    bestScore = -inf;

    % Build conflict map (logical sparse matrix)
    C = sparse(n, n);
    idx = sub2ind([n, n], min(CL, [], 2), max(CL, [], 2));
    C(idx) = true;
    C = C + C';  % make symmetric

    % Precompute preference list for each point
    [~, prefOrder] = sort(Xbar_shr, 2, 'descend');

    if verbose
        fprintf('Run\tObjective Value\t\tStatus\n');
    end
    for run = 1:numRuns
        X = false(n, k);
        clusterMembers = cell(k, 1);
        clusterCount = zeros(1, k);

        % Sort items by highest preference and shuffle
        [~, topChoice] = max(Xbar_shr, [], 2);
        order = randperm(n);
        [~, sortIdx] = sort(topChoice(order), 'descend');
        order = order(sortIdx);

        feasible = true;

        % Greedy assignment
        for pos = 1:n
            i = order(pos);
            assigned = false;
            for h = prefOrder(i, :)
                members = clusterMembers{h};
                if isempty(members)
                    % Always assign if cluster empty
                    X(i,h) = true;
                    clusterMembers{h}(end+1) = i;
                    clusterCount(h) = clusterCount(h) + 1;
                    assigned = true;
                    break;
                end
                % Check conflicts with assigned members
                if ~any(C(i, members))
                    X(i,h) = true;
                    clusterMembers{h}(end+1) = i;
                    clusterCount(h) = clusterCount(h) + 1;
                    assigned = true;
                    break;
                end
            end
            if ~assigned
                feasible = false;
                break;
            end
        end

        % Repair empty clusters
        if feasible
            emptyClusters = find(clusterCount == 0);
            for h = emptyClusters
                moved = false;
                for i = 1:n
                    currentH = find(X(i,:));
                    if isempty(currentH) || currentH == h
                        continue;
                    end
                    if ~any(C(i, clusterMembers{h}))
                        X(i,currentH) = false;
                        X(i,h) = true;
                        clusterMembers{currentH}(clusterMembers{currentH} == i) = [];
                        clusterMembers{h}(end+1) = i;
                        clusterCount(h) = 1;
                        clusterCount(currentH) = clusterCount(currentH) - 1;
                        moved = true;
                        break;
                    end
                end
                if ~moved
                    feasible = false;
                    break;
                end
            end
        end

        % Objective computation
        if feasible
            score = sum(Xbar_shr(X));
            if score > bestScore
                bestScore = score;
                X_best = X;
                if verbose
                    fprintf('%d\t%f\t\tfeasible\n', run, score);
                end
            end
        else
            bestScore = -inf;
            X_best = [];
            if verbose
                fprintf('%d\t--\t\tinfeasible\n', run);
            end
        end
    end

    tEnd = toc(tStart);

    X = X_best;
    f = bestScore;
    time = tEnd;
end