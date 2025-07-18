clc; clear; close all;

% [edit] directory to the project
path_root    = 'D:\Documents\BELANDA\PhD Thesis\Code\MATLAB\amode_navigation_experiment\experiment_b';

% declare some of the important paths
path_function = fullfile(path_root, 'functions');
path_outputs  = fullfile(path_root, 'outputs');
path_data     = fullfile(path_root, 'data');

addpath(genpath(path_function));

% [edit] directory to the trial
dir_trial    = "trial_0021_Session4_03";

% declare again some of the important paths
path_bonescan     = fullfile(path_data, dir_trial, 'bonescan');
path_intermediate = fullfile(path_data, dir_trial, 'intermediate');
path_measurement  = fullfile(path_data, dir_trial, 'measurement');

% get all the folders of the intermediates
items_dir = dir(path_intermediate);
folders_intermediate = {items_dir([items_dir.isdir] & ~ismember({items_dir.name}, {'.', '..'})).name};

% get all the folders of the measurements
items_dir = dir(path_measurement);
folders_measurement = {items_dir([items_dir.isdir] & ~ismember({items_dir.name}, {'.', '..'})).name};
% the last folder for measurements is always empty, we don't need that
folders_measurement(end) = [];


%%

% initialize 
ust_data = struct();

% get some constants
ust_dataspec.n_ust    = 30;
ust_dataspec.n_sample = 3500;

% loop for all folder measurements
for folder_idx = 1:length(folders_measurement)

    % get the current folder of intermediate
    folder_intermediate = fullfile(path_intermediate, folders_intermediate{folder_idx});
    % read the ultrasound data from intermediate
    [ust_datavalues_intermediate, ust_datatimestamps_intermediate, ust_datalabels_intermediate] = readTIFF_USsignal_v2(folder_intermediate, ust_dataspec.n_ust, ust_dataspec.n_sample, true);
    n_timestamp_intermediate = size(ust_datavalues_intermediate, 3);

    % get the current folder of measurement
    folder_measurement = fullfile(path_measurement, folders_measurement{folder_idx});
    % read the ultrasound data from measurement
    [ust_datavalues_measurement, ust_datatimestamps_measurement, ust_datalabels_measurement] = readTIFF_USsignal_v2(folder_measurement, ust_dataspec.n_ust, ust_dataspec.n_sample, false);
    n_timestamp_measurement = size(ust_datavalues_measurement, 3);

    % store the data
    ust_data.number = folder_idx-1;
    ust_data.intermediate.path       = folder_intermediate;
    ust_data.intermediate.values     = ust_datavalues_intermediate;
    ust_data.intermediate.timestamps = ust_datatimestamps_intermediate;
    ust_data.intermediate.labels     = ust_datalabels_intermediate;
    ust_data.measurement.path        = folder_measurement;
    ust_data.measurement.values      = ust_datavalues_measurement;
    ust_data.measurement.timestamps  = ust_datatimestamps_measurement;
    ust_data.measurement.labels      = ust_datalabels_measurement;

    % save the data
    filename = sprintf('ust_data_%s.mat',  folders_intermediate{folder_idx});
    save(fullfile(path_data, filename), 'ust_data', '-v7.3');

end