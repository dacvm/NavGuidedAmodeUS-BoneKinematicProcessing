%% HELOOOOOO
% - Before running this script, you must already generated .mat file called
%   "all_TsReg_sX_mXX.mat". This .mat file contains transformations of the
%   bone (tibia) relative to the ref (global) both GT and est.
% - The purpose of this script is to quantify and evaluate the 3d pose 
%   estimation against the ground truth.

clc; clear; close all;

% [EDIT] directory to the project
path_root    = 'D:\DennisChristie\NavGuidedAmodeUS-BoneKinematicProcessing';

% [EDIT] directory to the trial
dir_trial    = "trial_0023_Session4_04";

% [EDIT] Change the data you are using accordingly. 
% -----> dir_depthdata is created by main1_processDepthData.m
% -----> dir_Tdata is created by main3_registrationWithTime.m
dir_depthdata = 'depthdata_s4_m04_20250708-172830';
dir_Tdata     = 'Tdata_s4_m04_20250710-074800';

% [EDIT] Change the cycle timestamp name according to groundtruth from
% -----  which experiment you use.
cycle_timestamp_filename = 'cycle_timestamp_s4_m04.csv';

% [EDIT] for saving the resulting mat file
is_saveMat = true;


%% INITIALIZE PATHS AND LOADING SOME CONFIGURATION

% declare some of the important paths
path_function     = fullfile(path_root, 'functions');
path_outputs      = fullfile(path_root, 'outputs');
path_data         = fullfile(path_root, 'data');

% declare the path which consists of processed depth data
path_outputdepth = fullfile(path_outputs, 'output_allest', dir_depthdata);

% Generate path to function directory
addpath(genpath(path_function));



%% INITIALIZE SOME DATA

% 1) Load all the data related to bone pose estimation performance
% -> This mat file is generated from main3_registrationWithTime.m
tmp_str = split(dir_trial, '_');
sess_str = tmp_str{3};
meas_str = tmp_str{4};
mat_filename = sprintf('all_TsReg_s%s_m%s.mat', sess_str(end), meas_str);
mat_fullpath = fullfile(path_outputdepth, dir_Tdata, mat_filename);
load(mat_fullpath);

% 2) Load all the data related to knee joint 6 dof estimation.
% -> This mat file is generated from main4_kinematicEstimation.m
mat_filename = sprintf('all_kneeJoint6DOFs_s%s_m%s.mat', sess_str(end), meas_str);
mat_fullpath = fullfile(path_outputdepth, dir_Tdata, mat_filename);
load(mat_fullpath);

% 3) Load the cycle timestamp data. 
% -> This data indicates where the cycle motion starts, what was been 
%    done in the experiment.
% -> this is obtained from extra_detectFErotCycle.m
csv_filename = sprintf('cycle_timestamp_s%s_m%s.csv', sess_str(end), meas_str);
csv_fullpath = fullfile(path_outputs, csv_filename);
cycle_timestamp = readmatrix(csv_fullpath);

% 4) Alternative color scheme for some plots
c =  [116, 185, 255;
      255, 118, 117;
      253, 203, 110;
       85, 239, 196;
      253, 121, 168;
      162, 155, 254;
      178, 190, 195;
      129, 236, 236] / 255;  


%% INITIALIZE FIGURE OBJECTS AND EVERYTHING RELATED TO PLOTS

% First figure is to show the joint kinematic with all of the cycle parts
fig1 = figure('Name', 'Joint Kinematic: All Cycle Parts', 'Position', [100 100 1200 1000]);
ax1 = axes(fig1);
hold(ax1, 'on');
grid(ax1, 'on');
axis(ax1, 'equal');
axis(ax1, 'tight');
title(ax1, '3D');
xlabel(ax1, 'X');
ylabel(ax1, 'Y');
zlabel(ax1, 'Z');
view(ax1,[60, 10]);
ax1.FontSize = 12;


%% Figure 1

cycle_select     = [1,2];
cycle_idx_select = [cycle_timestamp(cycle_select(1)), cycle_timestamp(cycle_select(2))];

Ts_boneGT_ref        = all_TsReg_table.Ts_boneGT_ref( cycle_idx_select(1):cycle_idx_select(2) );
Ts_boneEst_ref       = all_TsReg_table.Ts_boneEstSmooth_ref( cycle_idx_select(1):cycle_idx_select(2) );
timestamp_idcs_valid = all_TsReg_table.Timestamp_idx( cycle_idx_select(1):cycle_idx_select(2) );
timestamp_ms_valid   = all_TsReg_table.Timestamp_ms( cycle_idx_select(1):cycle_idx_select(2) );
n_data_valid = length(timestamp_idcs_valid);

% convert cell array to 3d matrix, each slice is the T
Ts_boneGT_ref_mat = cat(3, Ts_boneGT_ref{:});
ts_boneGT_ref = squeeze(Ts_boneGT_ref_mat(1:3, 4, :));

% convert cell array to 3d matrix, each slice is the T
Ts_boneEstSmooth_ref_mat = cat(3, Ts_boneEst_ref{:});
ts_boneEstSmooth_ref = squeeze(Ts_boneEstSmooth_ref_mat(1:3, 4, :));

is_colorful = true;
if(is_colorful)

    % GT
    cmap = winter(n_data_valid);
    for i = 1:(n_data_valid-1)
        % plot segment [i  i+1]
        plot3( ts_boneGT_ref(1,i:i+1), ts_boneGT_ref(2,i:i+1), ts_boneGT_ref(3,i:i+1), '-', ...
               'Color'           , cmap(i,:), ...
               'MarkerFaceColor' , cmap(i,:), ...
               'MarkerEdgeColor' , cmap(i,:) );
    end
    
    idx = 1:5:n_data_valid-1;  
    scatter3( ts_boneGT_ref(1,idx), ts_boneGT_ref(2,idx), ts_boneGT_ref(3,idx), ...
              20,  cmap(idx,:), ...
              'o', 'filled' );

    % ESt
    cmap = autumn(n_data_valid);
    for i = 1:(n_data_valid-1)
        % plot segment [i  i+1]
        plot3( ts_boneEstSmooth_ref(1,i:i+1), ts_boneEstSmooth_ref(2,i:i+1), ts_boneEstSmooth_ref(3,i:i+1), '-', ...
               'Color'           , cmap(i,:), ...
               'MarkerFaceColor' , cmap(i,:), ...
               'MarkerEdgeColor' , cmap(i,:) );
    end
    
    idx = 1:5:n_data_valid-1;  
    scatter3( ts_boneEstSmooth_ref(1,idx), ts_boneEstSmooth_ref(2,idx), ts_boneEstSmooth_ref(3,idx), ...
              20,  cmap(idx,:), ...
              'o', 'filled' );

else

    % GT
    plot3(ts_boneGT_ref(1,:), ts_boneGT_ref(2,:), ts_boneGT_ref(3,:), '-b'); grid on; axis equal; hold on;
    idx = 1:5:n_data_valid;  
    scatter3( ts_boneGT_ref(1,idx), ts_boneGT_ref(2,idx), ts_boneGT_ref(3,idx), ...
              20, 'ob', 'filled' );

    % Est
    plot3(ts_boneEstSmooth_ref(1,:), ts_boneEstSmooth_ref(2,:), ts_boneEstSmooth_ref(3,:), '-r'); grid on; axis equal; hold on;
    idx = 1:5:n_data_valid;  
    scatter3( ts_boneEstSmooth_ref(1,idx), ts_boneEstSmooth_ref(2,idx), ts_boneEstSmooth_ref(3,idx), ...
              20, 'or', 'filled' );
end