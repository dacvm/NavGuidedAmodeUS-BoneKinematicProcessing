function estimatedbone_proc = mmodeDepthDetection_v2(mmode, depth_unit, seed, peakparam, movwinparam, clusterparam, outlierparam, smoothingparam, axis_result, display_process)
% MMODEDEPTHDETECTION  Detect bone surface depth in an M-mode ultrasound image
%
%   estimatedbone_proc = mmodeDepthDetection(mmode, depth_unit, new_seed, ...
%                                           peakparam, movwinparam, ...
%                                           clusterparam, outlierparam, ...
%                                           smoothingparam, axis_result, ...
%                                           display_process)
%
%   Summary
%   -------  
%   Processes a column-wise M-mode ultrasound image in chunks (moving window),
%   finds intensity peaks, clusters them into a bone surface estimate using
%   estimateBoneCluster_v3, and returns per-frame median depth and spread.
%
%   Inputs
%   ------  
%   mmode            Processed M-mode image (rows = depth samples, columns = time frames)
%   depth_unit       Depth represented by each row (in mm)
%   seed             Initial depth guess (in mm) for the clustering algorithm
%   peakparam        Struct with fields:
%                      > peak_threshold: minimum peak height
%                      > prom_threhsold: minimum peak prominence
%   movwinparam      (optional) Struct with fields:
%                      > length_window: window width (samples), default 500
%                      > length_stride: step size between windows, default 100
%   clusterparam     Parameters passed to estimateBoneCluster_v3 for clustering
%   outlierparam     Parameters for outlier rejection inside estimateBoneCluster_v3
%   smoothingparam   Smoothing settings for the clustering routine
%   axis_result      (optional) Handle to axes showing the mmode image;
%                    if non-empty, plots each window’s fitted bone line
%   display_process  Boolean flag; if true, shows and then closes clustering figures
%
%   Output
%   ------  
%   estimatedbone_proc  [K×3] array, where K = number of time frames:
%                       > Column 1: frame index (timestamp)
%                       > Column 2: median detected depth (mm)
%                       > Column 3: standard deviation of detected depths (mm)
%
%   How it works
%   -------------  
%   1. Divide the full mmode matrix into overlapping windows (centered on the start,
%      advancing by length_stride until the end).
%   2. In each window, for every column:
%        > Extract the 1-D depth signal.
%        > Find peaks above peak_threshold and with given prominence.
%        > Convert their row indices to depth in mm.
%   3. Collect all peaks in that window and call estimateBoneCluster_v3 with the
%      previous seed to cluster and fit a bone-surface curve.
%   4. Record the fitted depths (fitted curve evaluated at each frame) and update
%      new_seed to the mean depth of the current window for the next iteration.
%   5. (Optional) Overlay each window’s fit on axis_result, and display intermediate
%      clustering figures if requested.
%   6. After sliding through all windows, sort all detected points by frame index,
%      compute per-frame median and standard deviation, and return them.
%   7. Finally, if axis_result is provided, draw a shaded error band around the
%      median depth curve using plotShadedError.
%
%   Usage example
%   -------------  
%     % Prepare parameters
%     peakparam = struct('peak_threshold', 0.5, 'prom_threhsold', 0.2);
%     movwinparam = struct('length_window', 400, 'length_stride', 80);
%     % Call the function (no plotting)
%     result = mmodeDepthDetection(mmodeData, 0.01, 20, peakparam, ...
%                                  movwinparam, clusterparam, ...
%                                  outlierparam, smoothingparam, [], false);
%
%   See also estimateBoneCluster_v3, plotShadedError, adaptiveThreshInit

% if new_seed is a struct (means that the user should give lb and ub)
if(isstruct(seed))
    if(~any(isfield(seed, {'lb', 'm', 'ub'})))
        error('Error occured in in mmodeDepthDetection function, new_seed should be a scalar or a struct with (lb, m, ub) fields');
    end
    new_seed = seed.m;
else
    new_seed = seed;
end

% if moving window parameters are not present, let's make the default value
if(isempty(movwinparam))
    movwinparam.length_window = 500;
    movwinparam.length_stride = 100;
end

% get the length data
length_data       = size(mmode, 2);
length_window     = movwinparam.length_window;
length_stride     = movwinparam.length_stride;

% Start position: Move back by half window to center first window on data start
% First window will be centered around index 1
length_halfwindow = length_window/2;
idx_startwindow   = 1 - length_halfwindow + 1;
idx_window        = 1;

% variable to store the resulting estimated bone depth data
estimatedbone_all = [];

% Start depth detection with moving window


% Continue until the window center reaches near the end of data
% Stop when start_idx > (length_data - length_halfwindow) to ensure 
% last window is centered around the end of the data
while idx_startwindow <= length_data - length_halfwindow

     % Calculate theoretical end index of current window
    idx_endwindow = idx_startwindow + length_window - 1;
    
    % Handle boundaries. Some windows will extend beyond data boundaries, 
    % so we need to clip them
    idx_actualstart = max(1, idx_startwindow);
    idx_actualend   = min(length_data, idx_endwindow);
    
    % variable to store all detected peaks
    allpeaks_mm     = [];

    % Extract the actual window data (only the portion within data bounds)
    for idx_timestamp = idx_actualstart:idx_actualend

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

    % If the user specified the seed as a struct (which has lb and ub) it
    % means that the user want to focus on that region of interest.
    if(isstruct(seed))
        tmp_idx = find( (allpeaks_mm(:,2)>seed.lb) & (allpeaks_mm(:,2)<seed.ub) );
        allpeaksnew_mm = allpeaks_mm(tmp_idx, :);
    else
        allpeaksnew_mm = allpeaks_mm;
    end

    % estimate bone cluster
    [f, data_cleaned, ~, fig2] = estimateBoneCluster_v3( allpeaksnew_mm, new_seed, clusterparam, outlierparam, smoothingparam, display_process);
    
    % get the points based from the line estimation
    detectedbone_X    = data_cleaned(1,1):data_cleaned(end, 1);
    detectedbone_Y_mm = feval(f, detectedbone_X);
    estimatedbone_all = [estimatedbone_all; detectedbone_X', detectedbone_Y_mm];

    % Move to next window
    idx_startwindow = idx_startwindow + length_stride;
    idx_window      = idx_window + 1;
    new_seed        = mean(detectedbone_Y_mm);


    % display figure
    if(~isempty(axis_result))
        plot(axis_result, detectedbone_X, detectedbone_Y_mm, '-b', 'LineWidth', 1.5);
        drawnow;
    end

    % if user set to display the process, close it every end of the loop
    if(display_process)
        close(fig2);
    end

% end loop
end


% Structuring the depth data 

% sort it based on timestamp
estimatedbone_all = sortrows(estimatedbone_all, 1);

% calculate mean and std based on each individual timestamp
[unique_vals, ~, group_idx] = unique(estimatedbone_all(:,1));
means = splitapply(@median, estimatedbone_all(:,2), group_idx);
stds = splitapply(@std, estimatedbone_all(:,2), group_idx);

% Now uniqVals(i), means(i), stds(i) correspond to each unique first‐column value.
estimatedbone_proc = [unique_vals means stds];
display_shadedError(axis_result, estimatedbone_proc(:,1), estimatedbone_proc(:,2), estimatedbone_proc(:,3), 'colors', [1 0 0]);

end

