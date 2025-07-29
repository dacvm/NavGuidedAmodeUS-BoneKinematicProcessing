%% HELOOOO
% Summary of this script

clc; clear; close all;

% [EDIT] directory to the project
path_root  = 'D:\Documents\BELANDA\PhD Thesis\Code\MATLAB\amode_navigation_experiment\experiment_b';

% [EDIT] directory to the trial
dir_trials = {'trial_0011_Session3_02', 'trial_0023_Session4_04'};

% [EDIT] select which trial you want to use
trial_idx  = 1;
dir_trial  = dir_trials{trial_idx};

% [EDIT] Specify folder index (intermediate/measurement number)
folder_idx = 2;

% [EDIT] 
is_savemat = true;

%% INITIALIZE PATHS AND LOADING SOME CONFIGURATION

% declare some of the important paths
path_function     = fullfile(path_root, 'functions');
path_outputs      = fullfile(path_root, 'outputs');
path_data         = fullfile(path_root, 'data');

% declare some of important path of the data
path_bonescan     = fullfile(path_data, dir_trial, 'bonescan');
path_intermediate = fullfile(path_data, dir_trial, 'intermediate');
path_measurement  = fullfile(path_data, dir_trial, 'measurement');
path_bonestl  = fullfile(path_root, "data", "ct", "bone");

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

%%

% 1) Load the rigid body data that is stored in a .csv File
% -- Use the dir function to find all files with the .csv extension in the
%    directory. 
% -- This .csv file is the recording of the rigid bodies from Qualisys. 
% -- Here, we assume that there is always only one csv file the directory
fileList = dir(fullfile(path_measurement, folders_measurement{folder_idx}, '*.csv'));

% Check if the .csv file was found
if ~isempty(fileList)
    % Assuming there's only one CSV file, get its full name
    fullFileName = fullfile(path_measurement, folders_measurement{folder_idx}, fileList(1).name);
    % Read the CSV file and get the rigid bodies data (Table object)
    all_rigidbodies_table = readCSV_qualisysRigidBodies(fullFileName);
else
    % Display a message if no CSV file was found
    disp('No CSV file found in the specified directory.');
    return;
end

% save the mat file
if(is_savemat)
    % Get the necessary information about which data that we will use. 
    % I am taking the advantage of the naming format of the dir_trial
    tmp_str = split(dir_trial, '_');
    % get the 'sX' string
    sess_str = tmp_str{3};
    sess_str = ['s', sess_str(end)];
    % get the 'mXX' string
    meas_str = tmp_str{4};
    meas_str = ['m', meas_str];
    % get the 'dXX' string
    data_str = sprintf('d%02d', folder_idx);
    % format the naming
    mat_filename = sprintf('all_rigidbodies_table_%s_%s_%s.mat', sess_str, meas_str, data_str);
    save(mat_filename, 'all_rigidbodies_table');
end