function [Ixx, Iyy, Ixy] = ridgeDetection(I, varargin)
% RIDGEDETECTION Computes Hessian matrix components for ridge detection
%
% Syntax:
%   [Ixx, Iyy, Ixy] = ridgeDetection(I)                    % Default Gaussian smoothing (sigma=1.0)
%   [Ixx, Iyy, Ixy] = ridgeDetection(I, sigma)             % Gaussian smoothing with specified sigma
%   [Ixx, Iyy, Ixy] = ridgeDetection(I, 'none')            % No smoothing
%   [Ixx, Iyy, Ixy] = ridgeDetection(I, 'gaussian', sigma) % Explicit method and sigma
%
% Inputs:
%   I           - Input 2D matrix (grayscale image or data)
%   sigma       - Standard deviation for Gaussian smoothing (default: 1.0)
%   'none'      - No smoothing applied
%   'gaussian'  - Gaussian smoothing method (default)
%
% Outputs:
%   Ixx - Second derivative in x-direction (d²I/dx²)
%   Iyy - Second derivative in y-direction (d²I/dy²)
%   Ixy - Mixed second derivative (d²I/dxdy)
%
% Examples:
%   I = imread('image.jpg'); I = double(rgb2gray(I));
%   [Ixx, Iyy, Ixy] = ridgeDetection(I);           % Default: Gaussian with sigma=1.0
%   [Ixx, Iyy, Ixy] = ridgeDetection(I, 2.5);      % Gaussian with sigma=2.5
%   [Ixx, Iyy, Ixy] = ridgeDetection(I, 'none');   % No smoothing

    % Input validation
    if nargin < 1
        error('Input matrix I is required');
    end
    
    % Parse input arguments
    sigma = 1.0;      % Default sigma
    method = 'gaussian'; % Default method
    
    if nargin == 2
        if ischar(varargin{1}) || isstring(varargin{1})
            % Second argument is method string
            method = char(varargin{1});
            if strcmp(method, 'gaussian')
                sigma = 1.0; % Use default sigma for Gaussian
            end
        elseif isnumeric(varargin{1})
            % Second argument is sigma (implies Gaussian method)
            sigma = varargin{1};
            method = 'gaussian';
        else
            error('Second argument must be a numeric sigma value or method string');
        end
    elseif nargin == 3
        % Both method and sigma specified
        if (ischar(varargin{1}) || isstring(varargin{1})) && isnumeric(varargin{2})
            method = char(varargin{1});
            sigma = varargin{2};
        else
            error('Usage: ridgeDetection(I, method, sigma) where method is string and sigma is numeric');
        end
    elseif nargin > 3
        error('Too many input arguments');
    end
    
    % Validate method
    if ~ismember(method, {'gaussian', 'none'})
        error('Method must be ''gaussian'' or ''none''');
    end
    
    % Validate sigma for Gaussian method
    if strcmp(method, 'gaussian') && (sigma <= 0)
        error('Sigma must be positive for Gaussian smoothing');
    end
    
    % Convert to double for precision
    I = double(I);
    
    % Apply smoothing if requested
    if strcmp(method, 'gaussian') && sigma > 0
        % Use MATLAB's built-in Gaussian filtering
        % imgaussfilt is optimized and handles edge cases better
        if exist('imgaussfilt', 'file') == 2
            % Image Processing Toolbox function (preferred)
            I = imgaussfilt(I, sigma);
        else
            % Alternative using fspecial + imfilter (older versions)
            h = fspecial('gaussian', 2*ceil(3*sigma)+1, sigma);
            I = imfilter(I, h, 'same', 'replicate');
        end
    end
    
    % Compute first derivatives using MATLAB's gradient function
    [Ix, Iy] = gradient(I);
    
    % Compute second derivatives by applying gradient again
    [Ixx, Ixy_from_Ix] = gradient(Ix);
    [Ixy_from_Iy, Iyy] = gradient(Iy);
    
    % Average the two estimates of Ixy for better accuracy
    % (since ∂²I/∂x∂y should equal ∂²I/∂y∂x by Schwarz's theorem)
    Ixy = 0.5 * (Ixy_from_Ix + Ixy_from_Iy);
end