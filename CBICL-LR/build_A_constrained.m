function [A, b, n_CL] = build_A_constrained(original_n, k, init_ML, init_CL)

    % transformation for shrinking
    [T, ML_graph] = build_T(original_n, init_ML);
    n = size(T, 1); % new size
    CL = update_CL(ML_graph, init_CL);
    n_CL = size(CL, 1);

    Tt = T';
    e_n = ones(n, 1);
    e_s = T*ones(original_n, 1);
    TTt = T*Tt;
    
    % number of constraints
    Acell = cell(1, n+1+n_CL);

    for i=1:n
        Ai = sparse([repelem(i, n), 1:n], [1:n, repelem(i, n)], [0.5 .* e_n .* e_s; 0.5 .* e_n .* e_s], n, n);
        Acell{i} = Ai;
    end

    Acell{n+1} = TTt;

    for c = 1:n_CL
        i = CL(c, 1);
        j = CL(c, 2);
        Aij = sparse([i, j], [j, i], [0.5, 0.5], n, n);
        Acell{n+1+c} = Aij;
    end

    %keyboard


    % Assume C is your cell array of m sparse matrices, each of size n x n
    m = numel(Acell);
    
    % First, compute total number of nonzeros across all matrices
    nnz_total = sum(cellfun(@nnz, Acell));
    
    % Preallocate index and value arrays for sparse triplets
    I = zeros(nnz_total, 1);  % row indices
    J = zeros(nnz_total, 1);  % column indices
    V = zeros(nnz_total, 1);  % values
    
    idx = 1;
    
    for i = 1:m
        [r, c, v] = find(Acell{i});  % get non-zero entries of matrix i
        j = sub2ind([n n], r, c);  % linear indices in column-major order
    
        k = numel(v);
        I(idx:idx+k-1) = i;      % row i in A
        J(idx:idx+k-1) = j;      % column in A (flattened index)
        V(idx:idx+k-1) = v;      % value
    
        idx = idx + k;
    end
    
    % Now build sparse matrix A
    A = sparse(I, J, V, m, n^2);

    b = [e_n; k; zeros(n_CL, 1)];


end