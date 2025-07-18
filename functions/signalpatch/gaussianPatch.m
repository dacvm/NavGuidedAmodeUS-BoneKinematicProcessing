function patch = gaussianPatch(A, patchCenter, patchSize, sigma)
% GAUSSIANPATCH Extract a Gaussian-weighted patch from a matrix
%
% Inputs:
%   A         - Input matrix
%   row       - Center row position
%   col       - Center column position
%   patchSize - Size of the patch (scalar for square, [height width] for rectangular)
%   sigma     - Standard deviation of Gaussian:
%               - scalar: isotropic Gaussian (same sigma for both axes)
%               - [sigmaY, sigmaX]: anisotropic Gaussian (different sigma for each axis)
%               - optional, default: patchSize/6
%
% Output:
%   patch     - Gaussian-weighted patch

    % Handle input arguments
    if nargin < 3
        if isscalar(patchSize)
            sigma = patchSize / 6;  % Default sigma as 1/6 of patch size
        else
            sigma = min(patchSize) / 6;
        end
    end
    
    % Handle sigma input (isotropic vs anisotropic)
    if isscalar(sigma)
        sigmaY = sigma;  % Y-axis (rows) standard deviation
        sigmaX = sigma;  % X-axis (cols) standard deviation
    else
        sigmaY = sigma(1);  % Y-axis (rows) standard deviation
        sigmaX = sigma(2);  % X-axis (cols) standard deviation
    end
    
    % Handle patch size
    if isscalar(patchSize)
        patchHeight = patchSize;
        patchWidth = patchSize;
    else
        patchHeight = patchSize(1);
        patchWidth = patchSize(2);
    end
    
    % Get matrix dimensions
    [M, N] = size(A);
    
    % Calculate patch boundaries
    halfH = floor(patchHeight / 2);
    halfW = floor(patchWidth / 2);
    
    % Define patch region with boundary checking
    rowStart = max(1, patchCenter(1) - halfH);
    rowEnd = min(M, patchCenter(1) + halfH);
    colStart = max(1, patchCenter(2) - halfW);
    colEnd = min(N, patchCenter(2) + halfW);
    
    % Extract the raw patch
    rawPatch = A(rowStart:rowEnd, colStart:colEnd);
    
    % Create coordinate grids for the actual patch size
    actualHeight = rowEnd - rowStart + 1;
    actualWidth = colEnd - colStart + 1;
    
    % Create Gaussian weight matrix
    [X, Y] = meshgrid(1:actualWidth, 1:actualHeight);
    
    % Calculate center of the actual patch
    centerY = (actualHeight + 1) / 2;
    centerX = (actualWidth + 1) / 2;
    
    % Create anisotropic Gaussian weights
    gaussian = exp(-((X - centerX).^2 / (2 * sigmaX^2) + (Y - centerY).^2 / (2 * sigmaY^2)));
    
    % Apply Gaussian weighting
    patch = rawPatch .* gaussian;
end