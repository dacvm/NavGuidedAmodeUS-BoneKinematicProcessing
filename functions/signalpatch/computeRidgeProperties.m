function [ridgeStrength, orientationVector] = computeRidgeProperties(Ixx, Iyy, Ixy)
% COMPUTERIDGEPROPERTIES Computes ridge strength and orientation from Hessian
%
% Inputs:
%   Ixx, Iyy, Ixy - Hessian matrix components
%
% Outputs:
%   ridgeStrength     - Ridge strength (larger eigenvalue magnitude)
%   ridgeOrientation  - Ridge orientation in radians

    % Compute eigenvalues of Hessian matrix at each pixel
    % For 2x2 matrix [Ixx Ixy; Ixy Iyy], eigenvalues are:
    % lambda = 0.5 * (Ixx + Iyy ± sqrt((Ixx - Iyy)^2 + 4*Ixy^2))
    
    trace_H = Ixx + Iyy;
    det_H = Ixx .* Iyy - Ixy.^2;
    discriminant = sqrt((Ixx - Iyy).^2 + 4 * Ixy.^2);
    
    lambda1 = 0.5 * (trace_H + discriminant);
    lambda2 = 0.5 * (trace_H - discriminant);
    
    % Ridge strength is the magnitude of the larger eigenvalue
    ridgeStrength = max(abs(lambda1), abs(lambda2));
    
    % Ridge orientation (perpendicular to gradient direction)
    % ridgeOrientation = 0.5 * atan2(2 * Ixy, Ixx - Iyy);

    % Compute eigenvector components directly
    % For the more negative eigenvalue (lambda2), the eigenvector is:
    % If |Ixx - lambda2| > |Iyy - lambda2|, use first row: vx = Ixy, vy = lambda2 - Ixx
    % Otherwise use second row: vx = lambda2 - Iyy, vy = Ixy
    
    [rows, cols] = size(Ixx);
    orientationVector = zeros(rows, cols, 2);
    
    % Compute unnormalized eigenvector components
    diff_xx = abs(Ixx - lambda2);
    diff_yy = abs(Iyy - lambda2);
    
    % Choose more stable computation based on which difference is larger
    use_first_row = diff_xx >= diff_yy;
    
    % Method 1: Use first row of (H - λI)v = 0  =>  (Ixx-λ)vx + Ixy*vy = 0
    orientationVector(:,:,1) = use_first_row .* Ixy + (~use_first_row) .* (lambda2 - Iyy);
    orientationVector(:,:,2) = use_first_row .* (lambda2 - Ixx) + (~use_first_row) .* Ixy;
    
    % Normalize to unit vectors
    magnitude = sqrt(orientationVector(:,:,1).^2 + orientationVector(:,:,2).^2);
    
    % Avoid division by zero
    magnitude(magnitude < eps) = 1;
    
    orientationVector(:,:,1) = orientationVector(:,:,1) ./ magnitude;
    orientationVector(:,:,2) = orientationVector(:,:,2) ./ magnitude;
    
    % Note: These are unit vectors pointing in ridge direction
    % To verify: sqrt(orientationVector(:,:,1).^2 + orientationVector(:,:,2).^2) = 1
end