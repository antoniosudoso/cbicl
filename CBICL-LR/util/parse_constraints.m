function [ML_U, ML_V, CL_U, CL_V] = parse_constraints(filename)
    % Initialize arrays
    ML_U = [];
    ML_V = [];
    CL_U = [];
    CL_V = [];

    % Open file
    fid = fopen(filename, 'r');
    if fid == -1
        error('Cannot open the file.');
    end

    % Read each line
    while ~feof(fid)
        line = strtrim(fgetl(fid));
        if isempty(line)
            continue;
        end

        parts = strsplit(line);
        label = parts{1};
        u = str2double(parts{2});
        v = str2double(parts{3});

        % Append to corresponding arrays
        switch label
            case 'ML_U'
                ML_U(end+1, :) = [u, v];
            case 'ML_V'
                ML_V(end+1, :) = [u, v];
            case 'CL_U'
                CL_U(end+1, :) = [u, v];
            case 'CL_V'
                CL_V(end+1, :) = [u, v];
            otherwise
                warning('Unknown label: %s', label);
        end
    end

    fclose(fid);
end
