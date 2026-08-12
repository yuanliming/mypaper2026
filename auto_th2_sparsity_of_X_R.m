%% determine sparse structures of R and X according to theorem 2, and validate the results. 
clear; clc;

% 1. Define controller gain sparse structure mask (SK)
% SK is of size m x n, where m is input dimension, n is state dimension
SK = [1 1 0;
      0 1 1]; % this is input

[m, n] = size(SK);

% 2. Calculate the sparse mask for matrix X (mask_X)
% According to the zero-structure condition: X_ij = 0 if SK(k,i)=1 and SK(k,j)=0
SM = cell(1, m);
mask_X = ones(n, n); % Initialize mask_X as a full ones matrix (n x n)

for i = 1:m
    % Outer product of i-th row's non-zero pattern and zero pattern
    % Identifies prohibited off-diagonal entries in X
    SM{i} = SK(i, :)' * ~SK(i, :);
    
    % Element-wise update to set forbidden locations to 0
    mask_X = mask_X .* ~SM{i};
end

% The sparse structure mask for R directly matches SK
mask_R = SK;

% Display calculated mask matrices
disp('Calculated Sparse Mask for X (mask_X):');
disp(mask_X);

disp('Sparse Mask for R (mask_R):');
disp(mask_R);

% 3. Generate symbolic matrices R and X based on calculated masks
R = create_sparse_sym('R', mask_R);
X = create_sparse_sym('X', mask_X);

% --- Display Symbolic Matrices ---
disp('Symbolic Matrix R:');
disp(R);

disp('Symbolic Matrix X:');
disp(X);

% 4. Compute symbolic controller gain matrix K = R * inv(X)
K = R * inv(X);
disp('Resulting Controller Gain K (R * X^-1):');
disp(K);

%% Helper Function: Generate Symbolic Matrix from a 0-1 Mask Matrix
function S = create_sparse_sym(prefix, mask)
    % Get dimensions of the mask
    [rows, cols] = size(mask);
    
    % Initialize a symbolic zero matrix
    S = sym(zeros(rows, cols)); 
    
    % Loop through the mask and populate symbolic variables where mask == 1
    for i = 1:rows
        for j = 1:cols
            if mask(i, j) == 1
                % Create dynamic symbolic variable names (e.g., 'R_1_2', 'X_3_3')
                % Option: Add 'real' flag to enforce real-valued entries
                var_name = sprintf('%s_%d_%d', prefix, i, j);
                S(i, j) = sym(var_name, 'real');
            end
        end
    end
end