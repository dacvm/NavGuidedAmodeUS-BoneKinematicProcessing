function meanError = calculateRMSE(PC1, PC2, method)
% calculateMeanError Computes the mean distance error between two point clouds or meshes.
%
%   meanError = computeMeanError(PC1, PC2, method)
%
%   Inputs:
%       PC1    - Either an Nx3 matrix containing 3D points or a triangulation object representing a mesh.
%       PC2    - Triangulation object representing a complete mesh.
%       method - String, either 'pt2pt' for point-to-point or 'pt2plane' for point-to-plane error.
%
%   Output:
%       meanError - The computed mean error.
%
%   Example:
%       % If PC1 is an Nx3 matrix:
%       errorValue = computeMeanError(PC1, PC2, 'pt2plane');
%
%       % If PC1 is also a triangulation object:
%       errorValue = computeMeanError(PC1_mesh, PC2, 'pt2pt');

    % Check that PC2 is a triangulation object.
    if ~isa(PC2, 'triangulation')
        error('PC2 must be a triangulation object.');
    end

    % Determine if PC1 is a matrix or a triangulation object
    if isa(PC1, 'triangulation')
        points1 = PC1.Points;
    elseif ismatrix(PC1) && size(PC1, 2) == 3
        points1 = PC1;
    else
        error('PC1 must be either an Nx3 matrix or a triangulation object.');
    end

    % Validate the method argument.
    if nargin < 3
        error('You must provide PC1, PC2, and the method.');
    end
    if ~ismember(method, {'pt2pt', 'pt2plane'})
        error('Method must be either "pt2pt" or "pt2plane".');
    end

    % Find, for each point in points1, the nearest point in PC2.Points using knnsearch.
    [idx, dists] = knnsearch(PC2.Points, points1);
    
    if strcmp(method, 'pt2pt')
        % Point-to-point error: simply the Euclidean distance.
        meanError = mean(dists);
    else
        % Point-to-plane error: first compute vertex normals for PC2.
        normals = vertexNormal(PC2);
        
        % Compute difference vectors between points1 and the corresponding nearest points on PC2.
        diffVec = points1 - PC2.Points(idx, :);
        
        % Compute the component of each difference vector along the corresponding normal.
        pt2planeErrors = abs(sum(diffVec .* normals(idx, :), 2));
        
        % Mean point-to-plane error.
        meanError = mean(pt2planeErrors);
    end
end
