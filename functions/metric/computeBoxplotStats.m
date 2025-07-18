function [dataMin, lowerWhisker, Q1, Q2, Q3, upperWhisker, dataMax] = computeBoxplotStats(data)
% computeBoxplotStats calculates boxplot statistics from a dataset.
% Accepts either a vector (one distribution) or a matrix (each column is one distribution).
%
% Returns:
%   - dataMin: Minimum value for each distribution (row vector)
%   - lowerWhisker: Smallest data point >= (Q1 - 1.5*IQR) for each distribution
%   - Q1: 25th percentile (first quartile) for each distribution
%   - Q2: Median (50th percentile) for each distribution
%   - Q3: 75th percentile (third quartile) for each distribution
%   - upperWhisker: Largest data point <= (Q3 + 1.5*IQR) for each distribution
%   - dataMax: Maximum value for each distribution (row vector)
%
% Example usage:
%   data = randn(100,3); % Three distributions in 3 columns
%   [dataMin, lowerWhisker, Q1, Q2, Q3, upperWhisker, dataMax] = computeBoxplotStats(data);

    if isempty(data)
        error('Input data cannot be empty.');
    end

    % Determine the size: each column is a separate distribution.
    [n, m] = size(data);

    % Sort the data along each column.
    sortedData = sort(data, 1);

    % Compute min and max for each column.
    dataMin = min(sortedData, [], 1);
    dataMax = max(sortedData, [], 1);

    % Compute quartiles using prctile along the first dimension.
    Q1 = prctile(sortedData, 25, 1);
    Q2 = prctile(sortedData, 50, 1);
    Q3 = prctile(sortedData, 75, 1);

    % Compute the interquartile range (IQR) for each column.
    IQR = Q3 - Q1;

    % Compute potential whisker limits.
    lowerFence = Q1 - 1.5 * IQR;
    upperFence = Q3 + 1.5 * IQR;

    % Initialize whiskers.
    lowerWhisker = zeros(1, m);
    upperWhisker = zeros(1, m);

    % Determine actual whiskers for each distribution.
    for j = 1:m
        % For the lower whisker, find the smallest value in column j that is >= lowerFence.
        validLow = sortedData(:, j) >= lowerFence(j);
        lowerWhisker(j) = min(sortedData(validLow, j));

        % For the upper whisker, find the largest value in column j that is <= upperFence.
        validHigh = sortedData(:, j) <= upperFence(j);
        upperWhisker(j) = max(sortedData(validHigh, j));
    end
end
