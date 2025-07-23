function estimatedbone_proc = mmodeDepthDetection_v0(mmode, depth_unit, peakparam, clusterparam, outlierparam, smoothingparam, axis_result, display_process)

% get the length data
length_data       = size(mmode, 2);

% variable to store all detected peaks
allpeaks_mm     = [];

% loop over all timestamp
for idx_timestamp = 1:length_data

    % get the current signal
    current_signal = mmode(:, idx_timestamp);

    % find peaks
    [peak_val, peak_loc_idx]  = findpeaks(current_signal, "MinPeakHeight", peakparam.peak_threshold, 'MinPeakProminence', peakparam.prom_threhsold);

    % convert peak location from idx to mm (ust_usspec.ds)
    peak_loc_mm = peak_loc_idx * depth_unit;

    % store the peaks
    tmp_peaks_mm   = [idx_timestamp*ones(length(peak_loc_idx), 1) peak_loc_mm peak_val];
    allpeaks_mm    = [allpeaks_mm; tmp_peaks_mm];

end


% estimate bone cluster
[f, data_cleaned, ~, fig2] = estimateBoneCluster_v3( allpeaks_mm, [], clusterparam, outlierparam, smoothingparam, display_process);

% get the points based from the line estimation
detectedbone_X    = data_cleaned(1,1):data_cleaned(end, 1);
detectedbone_Y_mm = feval(f, detectedbone_X);
estimatedbone_all = [detectedbone_X', detectedbone_Y_mm];

% display figure
if(~isempty(axis_result))
    plot(axis_result, detectedbone_X, detectedbone_Y_mm, '-b', 'LineWidth', 1.5);
    drawnow;
end

% to be consistent with mmodeDepthDetection_v2. In that function, i include
% std (because i introduced moving window with overlap, such that each
% timestamp can have multiple detected depth value). Here we don't have
% stds, so let's put zero for that
estimatedbone_proc = [estimatedbone_all, zeros(length(detectedbone_Y_mm), 1)];

end

