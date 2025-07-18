%% HELOOOO
% Before running this script, you must convert the raw files from files in
% the data directory to mat files. This is just a convenient way to
% organize the data. Refer to: showMmodeWholeTrial_generateMAT.m

clc; clear; close all;

% [EDIT] directory to the project
path_root    = 'D:\Documents\BELANDA\PhD Thesis\Code\MATLAB\amode_navigation_experiment\experiment_b';

% [EDIT] directory to the trial
dir_trial    = "trial_0023_Session4_04";

% [EDIT] window configuration file
csvfile_windowconfig = 'transducerconfig_v8a_window2024-12-20_14-37-59_edited2025-07-03_15-43-26.csv';

% [EDIT] holder configuration file
csvfile_holderconfig = 'transducerconfig_v8a.csv';

% [EDIT] Specify folder index
folder_idx = 1;

% [EDIT] Specify if you want to load the intermediate data too alongside
%        the measurement data
is_withIntermediate = false;

% [EDIT] 
is_saveMat = false;


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

% 6) Get the session name
tmp_str = split(dir_trial, '_');
sess_str = tmp_str{3};
meas_str = tmp_str{4};

mat_filename = sprintf('all_depthgt_s%s_m%s.mat', sess_str(end), meas_str);
mat_fullpath = fullfile(path_outputs, 'output_allgtdepths', mat_filename);
load(mat_fullpath);


%% INITIALIZE CONSTANTS FOR ULTRASOUND PROCESSING

% get some constants
ust_dataspec.n_ust    = length(ust_windowconfig);
ust_dataspec.n_sample = 3500;

% select probes that has window defined
ust_config_selected = ust_windowconfig([ust_windowconfig.IsSet] == 1);
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
    n_timestamp = length(ust_data.intermediate.timestamps) + length(ust_data.measurement.timestamps);
else
    n_timestamp = length(ust_data.measurement.timestamps);
end

% This variables will be temporary output which stores all values of depth
% from all ultrasound transducer
ust_depthest_matrix = zeros(n_timestamp, n_ust_selected);
ust_depthstd_matrix = zeros(n_timestamp, n_ust_selected);

% variable to store all the ust number (index) we used
ust_numbers = zeros(1, n_ust_selected);

% Process and display m-mode image for all selected ust
for ust_idx=1:n_ust_selected

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

    % Display the window
    lb = ust_config_selected(ust_idx).LowerBound;
    m  = ust_config_selected(ust_idx).Middle;
    ub = ust_config_selected(ust_idx).UpperBound;
    % display_mmode_window(ax1{ust_idx}, lb, m, ub);
    plot(ax1{ust_idx}, 0, m, 'or', 'MarkerFaceColor', 'r', 'MarkerSize',12);
    rectangle(ax1{ust_idx}, 'Position', [0, lb, 20, ub-lb], 'FaceColor', 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.5);


    % [5] -----------------------------------------------------------------
    % Display ground truth
    timestamp_idx    = all_depthgtmean_table.("Timestamp_idx");
    ust_gtdepth_mean = all_depthgtmean_table.(num2str(ust_number));
    ust_gtdepth_std  = all_depthgtstd_table.(num2str(ust_number));
    display_shadedError(ax1{ust_idx}, timestamp_idx, ust_gtdepth_mean, ust_gtdepth_std, 'colors', [1 0 0]);


% end of ust_idx loop    
end