clear; clc; close all;

matfile = 'all_rigidbodies_table_s3_m02_d01.mat';
load(matfile);

% specify which column that i want to monitor
colNames = {'A_T_LEP', 'A_T_MAL', 'A_T_MEP', 'A_T_MID'};
% copy to another table
tmp_table = all_rigidbodies_table(:, colNames);

% get the number of rows and column
n_rows = height(tmp_table);
n_cols = width(tmp_table);
hasNaN = false(n_rows, n_cols);               % pre‑allocate logical matrix

for i = 1:n_cols
    % Each entry is a 1×1 cell → inside is a struct with field .T (4×4)
    hasNaN(:,i) = cellfun(@(s) any(isnan(s.T(:))), ...
                          tmp_table{:,i});
end

% Convert logical → binary (0/1) if you prefer a numeric matrix
binaryMatrix = double(hasNaN);


% ------------------------------------------------------------------------
% This part is for showing the report if there is NaN of each column.
% This form is easier to understand

% pre‑allocate a cell array to hold the row indices
nanRows_perCol = cell(1, n_cols);

for col = 1:n_cols
    % find all rows in this column where there’s a 1 (i.e. a NaN)
    nanRows_perCol{col} = find(binaryMatrix(:,col));
end

% Example: to display
for col = 1:n_cols
    fprintf('Column %2d → rows with NaNs: %s\n', ...
            col, mat2str(nanRows_perCol{col}));
end

% This part is just for removing the row with NaN -------------------------

% logical mask of rows with any 1
rowsWithNaN_logic = any(binaryMatrix, 2);

% numeric array of row numbers
rowsWithNaN = find(rowsWithNaN_logic);

% delete the row with problematic values
all_rigidbodies_table(rowsWithNaN, :) = [];

% save
tmp_strs    = strsplit(matfile, '.');
new_matfile = [tmp_strs{1}, '_clean', '.', tmp_strs{2}];
save(new_matfile, 'all_rigidbodies_table');