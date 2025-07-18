function ridgeMask = nonMaximumSuppression(ridgeStrength, orientationVector, varargin)
% NONMAXIMUMSUPPRESSION Performs non-maximum suppression along ridge axis
%
% Syntax:
%   ridgeMask = nonMaximumSuppression(ridgeStrength, orientationVector)
%   ridgeMask = nonMaximumSuppression(ridgeStrength, orientationVector, threshold)
%   ridgeMask = nonMaximumSuppression(ridgeStrength, orientationVector, threshold, method)
%
% Inputs:
%   ridgeStrength     - Ridge strength matrix from computeRidgeProperties
%   orientationVector - Ridge orientation vectors (rows x cols x 2) from computeRidgeProperties
%   threshold         - Minimum ridge strength threshold (default: 0.1 * max(ridgeStrength(:)))
%   method            - Interpolation method: 'bilinear' (default) or 'nearest'
%
% Outputs:
%   ridgeMask - Binary mask where 1 indicates ridge pixels (local maxima along ridge direction)
%
% Description:
%   This function performs non-maximum suppression by checking if each pixel's
%   ridge strength is a local maximum along the ridge direction. The ridge
%   direction is perpendicular to the gradient direction and is given by the
%   orientation vector from computeRidgeProperties.
%
% Example:
%   I = imread('fingerprint.jpg'); I = double(rgb2gray(I));
%   [Ixx, Iyy, Ixy] = ridgeDetection(I);
%   [ridgeStrength, orientationVector] = computeRidgeProperties(Ixx, Iyy, Ixy);
%   ridgeMask = nonMaximumSuppression(ridgeStrength, orientationVector);
%   figure; imshow(ridgeMask); title('Ridge Detection Result');

    % Input validation
    if nargin < 2
        error('ridgeStrength and orientationVector are required');
    end
    
    [rows, cols] = size(ridgeStrength);
    if size(orientationVector, 1) ~= rows || size(orientationVector, 2) ~= cols || size(orientationVector, 3) ~= 2
        error('orientationVector must have size [rows, cols, 2] matching ridgeStrength');
    end
    
    % Parse optional arguments
    threshold = 0.1 * max(ridgeStrength(:)); % Default threshold
    method = 'bilinear'; % Default interpolation method
    
    if nargin >= 3 && ~isempty(varargin{1})
        threshold = varargin{1};
        if ~isscalar(threshold) || threshold < 0
            error('threshold must be a non-negative scalar');
        end
    end
    
    if nargin >= 4 && ~isempty(varargin{2})
        method = varargin{2};
        if ~ismember(method, {'bilinear', 'nearest'})
            error('method must be ''bilinear'' or ''nearest''');
        end
    end
    
    % Initialize output mask
    ridgeMask = false(rows, cols);
    
    % Apply threshold first to reduce computation
    candidatePixels = ridgeStrength > threshold;
    
    % Create coordinate grids
    [X, Y] = meshgrid(1:cols, 1:rows);
    
    % Process each candidate pixel
    for i = 1:rows
        for j = 1:cols
            if ~candidatePixels(i, j)
                continue;
            end
            
            % Get ridge direction (orientation vector)
            vx = orientationVector(i, j, 1);
            vy = orientationVector(i, j, 2);
            
            % Skip if orientation vector is too small (numerical issues)
            if sqrt(vx^2 + vy^2) < eps
                continue;
            end
            
            % Current pixel strength
            currentStrength = ridgeStrength(i, j);
            
            % Sample points along ridge direction (both sides)
            % Use unit step size in the direction of the orientation vector
            step_size = 1.0;
            
            % Forward direction
            x_forward = j + step_size * vx;
            y_forward = i + step_size * vy;
            
            % Backward direction  
            x_backward = j - step_size * vx;
            y_backward = i - step_size * vy;
            
            % Check bounds and interpolate ridge strength values
            strength_forward = 0;
            strength_backward = 0;
            
            % Forward sample
            if x_forward >= 1 && x_forward <= cols && y_forward >= 1 && y_forward <= rows
                if strcmp(method, 'bilinear')
                    strength_forward = interp2(X, Y, ridgeStrength, x_forward, y_forward, 'linear', 0);
                else % nearest
                    x_round = round(x_forward);
                    y_round = round(y_forward);
                    if x_round >= 1 && x_round <= cols && y_round >= 1 && y_round <= rows
                        strength_forward = ridgeStrength(y_round, x_round);
                    end
                end
            end
            
            % Backward sample
            if x_backward >= 1 && x_backward <= cols && y_backward >= 1 && y_backward <= rows
                if strcmp(method, 'bilinear')
                    strength_backward = interp2(X, Y, ridgeStrength, x_backward, y_backward, 'linear', 0);
                else % nearest
                    x_round = round(x_backward);
                    y_round = round(y_backward);
                    if x_round >= 1 && x_round <= cols && y_round >= 1 && y_round <= rows
                        strength_backward = ridgeStrength(y_round, x_round);
                    end
                end
            end
            
            % Check if current pixel is local maximum along ridge direction
            if currentStrength >= strength_forward && currentStrength >= strength_backward
                % Additional check: ensure it's a significant local maximum
                % (not just due to numerical noise)
                if currentStrength > max(strength_forward, strength_backward) + eps || ...
                   (abs(currentStrength - max(strength_forward, strength_backward)) < eps && ...
                    currentStrength > threshold)
                    ridgeMask(i, j) = true;
                end
            end
        end
    end
    
    % Optional: Post-processing to remove isolated pixels (morphological opening)
    % Uncomment the following lines if you want to clean up the result
    % se = strel('disk', 1);
    % ridgeMask = imopen(ridgeMask, se);
end