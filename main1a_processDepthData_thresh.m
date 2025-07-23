%% HELOOOO
% Before running this script, you must convert the raw files from files in
% the data directory to mat files. This is just a convenient way to
% organize the data. Refer to: main0_showMmodeWholeTrial_generateMAT.m

clc; clear; close all;

% [EDIT] directory to the project
path_root    = 'D:\Documents\BELANDA\PhD Thesis\Code\MATLAB\amode_navigation_experiment\experiment_b';

% [EDIT] directory to the trial
dir_trial    = "trial_0011_Session3_02";

% [EDIT] window configuration file
csvfile_windowconfig = 'transducerconfig_v8a_window2024-12-19_11-29-46_edited2025-07-22_10-31-50.csv';

% [EDIT] holder configuration file
csvfile_holderconfig = 'transducerconfig_v8a.csv';

% [EDIT] Specify folder index
folder_idx = 1;

% [EDIT] Specify if you want to load the intermediate data too alongside
%        the measurement data
is_withIntermediate = false;

% [EDIT] 
is_saveMat_ustreqdata  = false;
is_saveMat_ustprocdata = false;
is_saveMat_final       = true;

% [EDIT]
is_loadCyclicData = false;


%% INITIALIZE PATHS AND LOADING SOME CONFIGURATION

% declare some of the important paths
path_function     = fullfile(path_root, 'functions');
path_outputs      = fullfile(path_root, 'outputs');
path_data         = fullfile(path_root, 'data');
% declare some of important path of the data
path_bonescan     = fullfile(path_data, dir_trial, 'bonescan');
path_intermediate = fullfile(path_data, dir_trial, 'intermediate');
path_measurement  = fullfile(path_data, dir_trial, 'measurement');

% 1) Generate path to function directory
addpath(genpath(path_function));

% 2) Get all the folders of the intermediates
items_dir = dir(path_intermediate);
folders_intermediate = {items_dir([items_dir.isdir] & ~ismember({items_dir.name}, {'.', '..'})).name};

% 3) Get all the folders of the measurements
items_dir = dir(path_measurement);
folders_measurement = {items_dir([items_dir.isdir] & ~ismember({items_dir.name}, {'.', '..'})).name};
% the last folder for measurements is always empty (it was like that from
% the navigation system) we don't need that.
folders_measurement(end) = [];

% 4) Load window configuration file
csv_fullpath = fullfile(path_intermediate, csvfile_windowconfig);
ust_windowconfig = table2struct(readtable(csv_fullpath));

% 5) Load holder configuration file
csv_fullpath = fullfile(path_data, 'configs', csvfile_holderconfig);
ust_holderconfig = table2struct(readtable(csv_fullpath));
% Find indices of structs where Group is not equal to 0, means that we are
% not using the ust
indices = arrayfun(@(x) x.Group ~= 0, ust_holderconfig);
% Create a new array of structs excluding those with Group equal to 0
ust_holderconfig = ust_holderconfig(indices);


% 6) Create folder for output if user want to save
if (is_saveMat_ustreqdata || is_saveMat_ustprocdata || is_saveMat_final)
    % get the session name
    tmp_str = split(dir_trial, '_');
    sess_str = tmp_str{3};
    meas_str = tmp_str{4};

    % get the current time
    foldername_timestamp  = datestr(now, 'yyyymmdd-HHMMSS');
    % create a unique name
    foldername_withprefix = sprintf('depthdata_s%s_m%s_%s', sess_str(end), meas_str, foldername_timestamp);
    % Full path and creation
    path_outputdepth = fullfile(path_outputs, 'output_allest', foldername_withprefix);
    mkdir(path_outputdepth);
end


% Load the cycle timestamp data. 
% -> This data indicates where the cycle motion starts, what was been 
%    done in the experiment.
% -> this is obtained from extra_detectFErotCycle.m
if(is_loadCyclicData)
    % get the session name
    tmp_str = split(dir_trial, '_');
    sess_str = tmp_str{3};
    meas_str = tmp_str{4};
    csv_filename = sprintf('cycle_timestamp_s%s_m%s.csv', sess_str(end), meas_str);
    csv_fullpath = fullfile(path_outputs, csv_filename);
    cycle_timestamp = readmatrix(csv_fullpath);
end


%% INITIALIZE CONSTANTS FOR ULTRASOUND PROCESSING

% get some constants
ust_dataspec.n_ust    = length(ust_windowconfig);
ust_dataspec.n_sample = 3500;

% Select probes that has window defined
% % The block below is the ideal one, selecting the ust based on the field
% % called .IsSet, which defined true/false when we do the experiment.
% ust_config_selected = ust_windowconfig([ust_windowconfig.IsSet] == 1);
% n_ust_selected      = length(ust_config_selected);
% But i will use this block, because in Session 3 (without navigation) we
% use full femur and tibia while in Session 4 (with navigation) we only use
% tibia. Devils in the detail, but we decided to only use tibia
ust_config_selected = ust_windowconfig(14:26);
n_ust_selected      = length(ust_config_selected);

% declare the ultrasound spesification
ust_usspec.vsound   = 1494;                   % [m/s]
ust_usspec.f        = 50 * 1e6;               % [Hz]
ust_usspec.dt       = 1/ust_usspec.f;         % [s]
ust_usspec.nt       = ust_dataspec.n_sample;  % [unit]
ust_usspec.t_vector = 0:ust_usspec.dt:(ust_usspec.dt*(ust_usspec.nt-1));     % [[s]]
ust_usspec.ds       = 0.5 * ust_usspec.vsound * ust_usspec.dt * 1e3;         % [mm]
ust_usspec.s_vector = 0.5 * ust_usspec.vsound * ust_usspec.t_vector * 1e3;   % [[mm]]
ust_usspec.pitch    = 8;                     % [mm]
ust_usspec.width    = 6;                     % [mm]
ust_usspec.kerf     = 2;                     % [mm]



%% INITIALIZE FIGURE OBJECTS

% prepare the figure object
fig1 = figure('Name', 'M-Mode Images');
% create a tab group
tg1 = uitabgroup(fig1);
% create tab and axis
ax1 = {};
for ust_idx=1:n_ust_selected
    % make a new tab
    tb1 = uitab(tg1, 'Title', sprintf('%s (%d)', ust_config_selected(ust_idx).GroupName, ust_config_selected(ust_idx).Number));
    % set the properties of the axes
    ax_tmp = axes(tb1);
    title(ax_tmp, 'M-mode Image');
    ylabel(ax_tmp,'Depth');
    xlabel(ax_tmp,'Timestamp');
    hold(ax_tmp, 'on');
    set(ax_tmp, 'YDir', 'reverse'); % Optionally set y-axis direction to normal
    % store the axes
    ax1{ust_idx} = ax_tmp;
end



%% MAIN PROGRAM

% Get the corresponding data
fprintf('Loading ust_data ...'); tic;
% Here i used folders_intermediate, it doesn't matter because in the
% navigation software, intermediate folders is always created along side 
% with measurement folders.
filename = sprintf('ust_data_%s.mat', folders_intermediate{folder_idx});
fullpath = fullfile(path_data, 'matfiles', dir_trial, filename);
load(fullpath);
fprintf(' (%.2fs)\n', toc);

% Get the number of timestamp, it depends whether the user wants to include
% the intermediate data or not
if (is_withIntermediate)
    % calculate the number of timestamps we use
    n_timestamp = length(ust_data.intermediate.timestamps) + length(ust_data.measurement.timestamps);

    % convert the timestamps into seconds and then store it
    timestamps  = [ust_data.intermediate.timestamps; ust_data.measurement.timestamps];   % in [ms]
    timestamps  = timestamps - timestamps(1);

else
    n_timestamp = length(ust_data.measurement.timestamps);
    timestamps  = ust_data.measurement.timestamps;
    timestamps  = timestamps - timestamps(1);
end

% This variables will be temporary output which stores all values of depth
% from all ultrasound transducer
ust_depthest_matrix = zeros(n_timestamp, n_ust_selected);
ust_depthstd_matrix = zeros(n_timestamp, n_ust_selected);
ust_depthamp_matrix = zeros(n_timestamp, n_ust_selected);

% variable to store all the ust number (index) we used
ust_numbers = zeros(1, n_ust_selected);

%%

% set these thresholds to the number of the n_ust_selected pwease.....
noise_constant = 4;
noise_weight   = [1 1 1 1 1 1 1 1 1 1 1 1 1];
tresh_allconstants = noise_constant * noise_weight;

% Process and display m-mode image for all selected ust
for ust_idx=1:n_ust_selected
% for ust_idx=1

    % [1] -----------------------------------------------------------------
    % Get the current ust number
    ust_number           = ust_config_selected(ust_idx).Number;
    ust_numbers(ust_idx) = ust_number;
    % Get the measurement matrix
    mmode_matrix_measurement = squeeze(ust_data.measurement.values(ust_number, :, :));
    
    % [2] -----------------------------------------------------------------
    % If you specified this to true, you will load the intermediate data too
    if (is_withIntermediate)
        % For other folders: display both intermediate and measurement matrices side-by-side
        mmode_matrix_intermediate = squeeze(ust_data.intermediate.values(ust_number, :, :));
        mmode_matrix = [mmode_matrix_intermediate, mmode_matrix_measurement];
    
    % If not, you will only load the measurement data
    else
        % For first folder: display only the measurement matrix
        mmode_matrix = mmode_matrix_measurement;
    end

    % [3] -----------------------------------------------------------------
    % Filter the matrix with sigmoid function
    cutoff_row       = round(1.25 / ust_usspec.ds);
    midpoint_ratio   = 0.5;
    steepness        = 2.0;
    mmode_matrix = applyPartialSigmoid(mmode_matrix, cutoff_row, steepness, midpoint_ratio);


    % [4] -----------------------------------------------------------------
    % Adjust the parameters required by display_mmode()
    mmode_x           = [1, n_timestamp];
    mmode_y           = [ust_usspec.s_vector(1), ust_usspec.s_vector(end)];
    mmode_imagethresh = 10000;
    mmode_proc = display_mmode(ax1{ust_idx}, mmode_matrix, mmode_imagethresh, mmode_x, mmode_y);
    xline(ax1{ust_idx}, n_timestamp, '-k', 'Linewidth', 2);

    % If you load a cyclic data, let's show it in a line
    if(is_loadCyclicData)
        for i=1:length(cycle_timestamp)
            xline(ax1{ust_idx}, cycle_timestamp(i), '-g', 'LineWidth', 1);
        end
    end

    % Display the window
    lb = ust_config_selected(ust_idx).LowerBound;
    m  = ust_config_selected(ust_idx).Middle;
    ub = ust_config_selected(ust_idx).UpperBound;
    % display_mmode_window(ax1{ust_idx}, lb, m, ub);
    plot(ax1{ust_idx}, 0, m, 'or', 'MarkerFaceColor', 'r', 'MarkerSize',12);
    rectangle(ax1{ust_idx}, 'Position', [0, lb, 20, ub-lb], 'FaceColor', 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.5);


    % [5] -----------------------------------------------------------------
    % Structure the data (and save if user want)
    ust_datapacket.number     = ust_number;
    ust_datapacket.mmode_raw  = mmode_matrix;
    ust_datapacket.mmode_proc = mmode_proc;
    ust_datapacket.mmode_x    = mmode_x;
    ust_datapacket.mmode_y    = mmode_y;
    ust_datapacket.init       = m;

    % Save the data if you want, for documentation in the future
    if(is_saveMat_ustreqdata)
        mat_filename = sprintf('ust_reqdata_%02d.mat', ust_number);
        mat_fullpath = fullfile(path_outputdepth, mat_filename);
        save(mat_fullpath, 'ust_datapacket', 'ust_usspec', 'ust_dataspec');
    end

    % [6] -----------------------------------------------------------------
    % Initialize the threshold. The threshold will be defined as a constant
    % time the noise of a particular ultrasound transducer.
    noise_enddepth_mm  = 1; % mm
    noise_enddepth_idx = round(noise_enddepth_mm/ust_usspec.ds);
    noise_patch        = ust_datapacket.mmode_proc(end-noise_enddepth_idx:end, :);
    noise_mean         = mean(noise_patch, 'all');
    peak_threshold     = tresh_allconstants(ust_idx) * noise_mean;
    prom_threshold     = peak_threshold;


    % [7] -----------------------------------------------------------------

    % parameters for depth
    depth_unit = ust_usspec.ds;

    % peak threshold parameters
    peakparam.peak_threshold = peak_threshold;
    peakparam.prom_threhsold = prom_threshold;

    % bone cluster parameters

    if(ust_idx==6)
        clusterparam.epsilon    = 0.09;
        clusterparam.minpts     = 100;
        clusterparam.cov        = [300 0; 0 1.0];
        outlierparam.gradthresh = 5;
        outlierparam.minpoint   = 15;
        smoothingparam          = 0.0005;
    else
        clusterparam.epsilon    = 0.09;
        clusterparam.minpts     = 100;
        clusterparam.cov        = [300 0; 0 1.0];
        outlierparam.gradthresh = 4.5;
        outlierparam.minpoint   = 20;
        smoothingparam          = 0.0003;
    end
    % detect the bone from mmode image
    estimatedbone_proc = mmodeDepthDetection_v0( ust_datapacket.mmode_proc, depth_unit, peakparam, ...
                                                  clusterparam, outlierparam, smoothingparam, ...
                                                  ax1{ust_idx}, true);

    % [8] -----------------------------------------------------------------
    % Get the amplitude of the processed datapoints
    estimatedbone_depthidx   = round(estimatedbone_proc(:,2) / ust_usspec.ds); % row
    estimatedbone_timeidx    = estimatedbone_proc(:,1);                        % col
    idx = sub2ind(size(ust_datapacket.mmode_raw), estimatedbone_depthidx, estimatedbone_timeidx);
    estimatedbone_depthamp   = ust_datapacket.mmode_raw(idx);

    % Structure the data
    ust_depthest.number    = ust_number;
    ust_depthest.timestamp = estimatedbone_proc(:,1);
    ust_depthest.mean      = estimatedbone_proc(:,2);
    ust_depthest.std       = estimatedbone_proc(:,3);
    ust_depthest.amp       = estimatedbone_depthamp;

    % save the data here 
    if(is_saveMat_ustprocdata)
        % save
        mat_filename = sprintf('ust_procdata_%02d.mat', ust_number);
        mat_fullpath = fullfile(path_outputdepth, mat_filename);
        save(mat_fullpath, 'ust_depthest');
    end

    % [8] -----------------------------------------------------------------
    ust_depthest_matrix(ust_depthest.timestamp, ust_idx) = ust_depthest.mean;
    ust_depthstd_matrix(ust_depthest.timestamp, ust_idx) = ust_depthest.std;
    ust_depthamp_matrix(ust_depthest.timestamp, ust_idx) = ust_depthest.amp;


% end of ust_idx loop    
end

% initialize a table
ust_cell     = arrayfun(@num2str, ust_numbers, 'UniformOutput', false);
table_header = [{'Timestamp_idx', 'Timestamp_ms'}, ust_cell];

% convert matrix to table
ust_timestampidx_vector = 1:n_timestamp;
all_depthestmean_table  = array2table([ust_timestampidx_vector', timestamps, ust_depthest_matrix], 'VariableNames', table_header);
all_deptheststd_table   = array2table([ust_timestampidx_vector', timestamps, ust_depthstd_matrix], 'VariableNames', table_header);
all_depthestamp_table   = array2table([ust_timestampidx_vector', timestamps, ust_depthamp_matrix], 'VariableNames', table_header);

if(is_saveMat_final)
    % get the session name
    tmp_str = split(dir_trial, '_');
    sess_str = tmp_str{3};
    meas_str = tmp_str{4};
    % save
    mat_filename = sprintf('all_depthest_s%s_m%s.mat', sess_str(end), meas_str);
    mat_fullpath = fullfile(path_outputdepth, mat_filename);
    save(mat_fullpath, 'all_depthestmean_table', 'all_deptheststd_table', 'all_depthestamp_table');
end