function [peak_threshold, prom_threshold] = adaptiveThreshInit(mmmode, init_depth_idx, init_window_idx, display)
% adaptiveThreshInit  Initialize intensity thresholds for an M-mode ultrasound image
%
%   [peak_threshold, prom_threshold] = adaptiveThreshInit(mmmode, ...
%       init_depth_idx, init_window_idx, display)
%
%   This function examines the signal strength within a specified depth
%   window of an M-mode image and computes two adaptive thresholds:
%     • peak_threshold    – the scaled maximum of a weighted median signal
%     • prom_threshold    – a fraction of the peak_threshold, for prominence
%
%   Inputs:
%     mmmode           – 2D matrix representing the processed M-mode image
%     init_depth_idx   – scalar index (row) for the center of the depth window
%     init_window_idx  – half-width (in rows) of the depth window
%     display          – boolean flag; if true, plot the windowed signals
%
%   Outputs:
%     peak_threshold   – 0.7 × max(weighted median signal in window)
%     prom_threshold   – 0.3 × peak_threshold
%
%   Example:
%     % use a ±10-row window around depth row 50, and show plots
%     [pth, prm] = adaptiveThreshInit(mmode_img, 50, 10, true);
%
%   Notes:
%     – A Gaussian window (Gausswin) weights the median signal before
%       finding its maximum.
%     – The prominence threshold is set lower for downstream peak
%       prominence calculations.
%
%   See also median, gausswin, max, plot

    % Determine the row-range for the initialization window
    init_window_idcs = init_depth_idx + [-1, 1] * init_window_idx;
    % Clamp to image top if underflow
    if init_window_idcs(1) < 1
        init_window_idcs(1) = 1;
    end



    % Extract the slice of the M-mode image within our window
    signalclipped_mmode = mmmode(init_window_idcs(1) : init_window_idcs(2), :);

    % Compute the median intensity across columns, yielding a depth profile
    signalclipped_mean = median(signalclipped_mmode, 2);

    % Create a Gaussian weighting vector to emphasize central depths
    %   - length matches number of rows in our window
    %   - shape tuned by parameter (here 1.2)
    signalclipped_weight = gausswin(length(signalclipped_mean), 1.2);

    % Apply weighting to the median signal
    signalclipped_weighted = signalclipped_mean .* signalclipped_weight;

    % Find the maximum of the weighted median signal
    signalclipped_max = max(signalclipped_weighted);



    % Initialize peak threshold to 70% of that maximum
    peak_threshold = 0.7 * signalclipped_max;

    % Set prominence threshold to 30% of the peak threshold
    prom_threshold = 0.3 * peak_threshold;


    % Check the threshold if it falls below noise level
    noise = median(mmmode(end-100:end, :), 'all');
    if (peak_threshold < noise+(0.25*noise) )
        peak_threshold = noise+(0.25*noise);
    end


    % If requested, display all intermediate signals for inspection
    if display
        figure;
        hold on;
        for idx=1:size(signalclipped_mmode,1)
            plot(signalclipped_mmode(:, idx), '--', 'LineWidth', 0.25);
        end
        plot(signalclipped_mean, '-r', 'LineWidth', 2);
        plot(signalclipped_weighted, '-m', 'LineWidth', 2);
        plot(signalclipped_weight, '--b', 'LineWidth',2);
        axis tight;
    end

end
