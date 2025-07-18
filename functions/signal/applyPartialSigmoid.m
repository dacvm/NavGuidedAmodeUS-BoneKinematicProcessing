function [weighted_matrix, weights] = applyPartialSigmoid(input_matrix, N, steepness, midpoint_ratio)
% APPLY_SIGMOID_ROW_WEIGHTING Apply sigmoid-based row weighting to a matrix
%
% INPUTS:
%   input_matrix    - Input matrix to be weighted (any size)
%   N               - Number of rows over which to apply sigmoid weighting
%   steepness       - Controls sharpness of sigmoid transition (higher = sharper)
%                     Typical range: 0.5 to 3.0 (default: 1.0)
%   midpoint_ratio  - Ratio determining sigmoid midpoint within first N rows
%                     0.5 means midpoint at N/2 (default: 0.5)
%                     Range: 0.1 to 0.9
%
% OUTPUTS:
%   weighted_matrix - Input matrix with sigmoid weighting applied to rows
%   weights         - Weight vector used (sigmoid for first N rows, 1 for rest)
%
% EXAMPLE:
%   % Create test matrix
%   test_matrix = rand(3500, 2000);
%   
%   % Apply sigmoid weighting over first 300 rows
%   [filtered_matrix, weight_vector] = apply_sigmoid_row_weighting(test_matrix, 300, 0.5, 0.5);
%   
%   % Plot weights to visualize
%   figure; plot(weight_vector(1:500)); title('Row Weights');

    % Input validation
    if nargin < 2
        error('At least 2 inputs required: input_matrix and N');
    end
    if nargin < 3 || isempty(steepness)
        steepness = 1.0;  % Default steepness
    end
    if nargin < 4 || isempty(midpoint_ratio)
        midpoint_ratio = 0.5;  % Default midpoint at N/2
    end
    
    % Validate inputs
    if ~ismatrix(input_matrix)
        error('input_matrix must be a 2D matrix');
    end
    if N <= 0 || N ~= round(N)
        error('N must be a positive integer');
    end
    if steepness <= 0
        error('steepness must be positive');
    end
    if midpoint_ratio <= 0 || midpoint_ratio >= 1
        error('midpoint_ratio must be between 0 and 1');
    end
    
    % Get matrix dimensions
    [num_rows, num_cols] = size(input_matrix);
    
    % Ensure N doesn't exceed matrix rows
    N = min(N, num_rows);
    
    % Create row indices for the transition region
    row_indices = (1:N)';
    
    % Calculate sigmoid parameters for perfect continuity
    midpoint = midpoint_ratio * N;
    
    % Create sigmoid that goes from near 0 to near 1 over the range
    % Use a symmetric range around midpoint to ensure smooth transition
    half_range = steepness * N / 4;  % Control the transition width
    
    % Calculate raw sigmoid values
    raw_sigmoid = 1 ./ (1 + exp(-(row_indices - midpoint) / half_range * 6));
    
    % Force perfect continuity: normalize so that sigmoid(1) ≈ 0 and sigmoid(N) = 1
    sigmoid_at_1 = 1 / (1 + exp(-(1 - midpoint) / half_range * 6));
    sigmoid_at_N = 1 / (1 + exp(-(N - midpoint) / half_range * 6));
    
    % Linear rescaling to ensure: weight(1) starts near 0, weight(N) = exactly 1
    sigmoid_weights = (raw_sigmoid - sigmoid_at_1) / (sigmoid_at_N - sigmoid_at_1);
    
    % Create full weight vector
    weights = ones(num_rows, 1);
    weights(1:N) = sigmoid_weights;
    
    % Apply weights to matrix (broadcasting automatically handles row-wise multiplication)
    weighted_matrix = input_matrix .* weights;
    
    % Optional: Display information about the weighting
    if nargout == 0  % If no output arguments, display info
        fprintf('Sigmoid weighting applied:\n');
        fprintf('  Matrix size: %d x %d\n', num_rows, num_cols);
        fprintf('  Transition rows: %d\n', N);
        fprintf('  Weight range: [%.6f, %.6f]\n', min(weights), max(weights));
        fprintf('  Weight at row N: %.6f (should be 1.000000)\n', weights(N));
        fprintf('  Steepness: %.3f\n', steepness);
        fprintf('  Midpoint ratio: %.3f (row %.1f)\n', midpoint_ratio, midpoint);
        fprintf('  Continuity check: weight(%d) = %.6f, weight(%d) = %.6f\n', ...
                N, weights(N), min(N+1, num_rows), weights(min(N+1, num_rows)));
    end
end