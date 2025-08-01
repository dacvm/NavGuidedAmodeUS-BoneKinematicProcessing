%% HELOOOOOO
% - Before running this script, you must already generated .mat file called
%   "all_TsReg_sX_mXX.mat". This .mat file contains transformations of the
%   bone (tibia) relative to the ref (global) both GT and est.
% - The purpose of this script is to quantify and evaluate the 3d pose 
%   estimation against the ground truth.

clc; clear; close all;

% [EDIT] directory to the project
path_root    = 'D:\Documents\BELANDA\PhD Thesis\Code\MATLAB\amode_navigation_experiment\experiment_b';

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
fig1 = figure('Name', 'Joint Kinematic: All Cycle Parts', 'Position', [50 50 1200 500]);
t1 = tiledlayout(fig1, 1, 2, ...
     'TileSpacing', 'compact', ...   % tighten spacing if you like
     'Padding',     'compact');      % remove outer margins
% Pre-allocate an array of axes handles
ax1 = gobjects(2,1);

% plot 1
ax1(1) = nexttile(t1, 1);
hold(ax1(1), 'on');
grid(ax1(1), 'on');
axis(ax1(1), 'tight');
title(ax1(1), 'Mean Registration Error Through Time');
ax1(1).FontSize = 12;
% left
yyaxis(ax1(1),'left')
xlabel(ax1(1), 'Timestamp');
ylabel(ax1(1), 'Angle (deg)');
% right
yyaxis(ax1(1),'right')
ylabel(ax1(1), 'Reg Error (for all cycle)');


% plot 2
ax1(2) = nexttile(t1, 2);
hold(ax1(2), 'on');
grid(ax1(2), 'on');
axis(ax1(2), 'tight');
title(ax1(2), 'Registration error');
xlabel(ax1(2), 'Distance Metric');
ylabel(ax1(2), 'RMSE (mm)');
ax1(2).FontSize = 12;


%% Figure 1


% Grab the valid index
timestamp_idcs_valid = all_TsReg_table.Timestamp_idx;

% prepare the variable
all_RMSE_value  = [];
all_RMSE_grpidx = [];
all_RMSE_grpstr = {};
xtlabels_str = {'us2estbone (p2p)', 'us2estbone (p2pl)', 'estbone2gtbone (p2p)', 'estbone2gtbone (p2pl)'};

% get the withnav value
RMSE_value_winav  = [all_TsReg_table.errors_rmse_us2regbone, all_TsReg_table.errors_rmse_regbone2gtbone];
RMSE_grpidx_winav = 1 .* ones(1, size(RMSE_value_winav, 1));

% this is dummy, for without nav value
RMSE_value_wonav  = zeros(size(RMSE_value_winav));
RMSE_grpidx_wonav = 2 .* ones(1, size(RMSE_value_wonav, 1));

% structure the data 
all_RMSE_value  = [RMSE_value_wonav; RMSE_value_winav];
all_RMSE_grpidx = [RMSE_grpidx_wonav, RMSE_grpidx_winav];
all_RMSE_grpstr = {'Without Navigation', 'With Navigation'};

% activate the corresponding axis (daboxplot does not have argument to 
% select particular axes)
axes(ax1(2));
% draw the bar plot
dabarplot( all_RMSE_value, ...
            'groups', all_RMSE_grpidx, ...
            'legend', all_RMSE_grpstr, ...
            'xtlabels', xtlabels_str, ...
            'errorbars', 'SD',...
            'errorhats', 0, ...
            'scatter', 2, ...
            'scatteralpha', 0.1, ...
            'baralpha', 0.6, ...
            'colors', c);


%%  Figure 2

% Grab the valid index
timestamp_idcs_valid = all_kneeJoint6DOFs_table.Timestamp_idx;

% Grab the knee joint 6dof and convert it into matrix
rmse_us2regbone     = all_TsReg_table.errors_rmse_us2regbone;
rmse_regbone2gtbone = all_TsReg_table.errors_rmse_regbone2gtbone;
kneeJoint6DOFs_gt   = cell2mat(all_kneeJoint6DOFs_table.kneeJoint6DOFs_gt);

% I pick 4, this one is for flexion-extension. How do i know it, because i
% am the one who structured it. Yes it is arbitrary. Who cares. Now go.
selectDoF_kneeJoint_gt = kneeJoint6DOFs_gt(:,4);

% before we do further processing, i want to chop the knee joint data
% into different cyclic motion parts. here, i initialize a cell to 
% store all the parts (parts can have different length)
rmse_us2regbone_cycleparts         = {};
rmse_regbone2gtbone_cycleparts         = {};
selectDoF_kneeJoint_gt_cycleparts  = {};

% get the cycle timestamp
for i=1:(length(cycle_timestamp)-1)
    % get the start and end index
    idx_start = cycle_timestamp(i);
    idx_end = cycle_timestamp(i+1);
    
    % before we go further let's check if the end index is bigger than
    % the valid timestamp we have
    if(idx_end>timestamp_idcs_valid(end))
        % if yes, let's not use this specific cycle part, stop it.
        break;
    end

    % get the current indices of the cycle part
    selectDoF_cyclepart_idcs  = find( (timestamp_idcs_valid>=idx_start) & (timestamp_idcs_valid<idx_end) );

    % get the current value of the cycle part
    rmse_us2regbone_cyclepart_values     = rmse_us2regbone(selectDoF_cyclepart_idcs);
    rmse_regbone2gtbone_cyclepart_values = rmse_regbone2gtbone(selectDoF_cyclepart_idcs);
    selectDoF_cyclepart_gt_values        = selectDoF_kneeJoint_gt(selectDoF_cyclepart_idcs);

    % store the current values of this cycle part
    rmse_us2regbone_cycleparts{i}        = rmse_us2regbone_cyclepart_values;
    rmse_regbone2gtbone_cycleparts{i}    = rmse_regbone2gtbone_cyclepart_values;
    selectDoF_kneeJoint_gt_cycleparts{i} = selectDoF_cyclepart_gt_values;
end



% get how many parts we have (both gt and est must be the same)
n_cycleparts = length(selectDoF_kneeJoint_gt_cycleparts);
% get the max length of the vector inside the cell
max_length = max( cellfun(@numel, selectDoF_kneeJoint_gt_cycleparts) );

% allocate memory to store the new stretched version of the part
rmse_us2regbone_cyclepartsnew        = zeros(max_length, n_cycleparts);
rmse_regbone2gtbone_cyclepartsnew    = zeros(max_length, n_cycleparts);
selectDoF_kneeJoint_gt_cyclepartsnew = zeros(max_length, n_cycleparts);

% we will loop for each cycle parts and we will stretch the shorter 
% parts to be the same length as the longest part.
for i=1:n_cycleparts

    % get the current cycle part
    rmse_us2regbone_cyclepart_values     = rmse_us2regbone_cycleparts{i};
    rmse_regbone2gtbone_cyclepart_values = rmse_regbone2gtbone_cycleparts{i};
    selectDoF_cyclepart_gt_values        = selectDoF_kneeJoint_gt_cycleparts{i};

    % stretch the value, here i am using imresize from image processing
    % tool box, because why not? they already have built-in
    % interpolation function inside
    rmse_us2regbone_cyclepart_valuesstretched     = imresize(rmse_us2regbone_cyclepart_values, [max_length 1], 'bicubic');
    rmse_regbone2gtbone_cyclepart_valuesstretched = imresize(rmse_regbone2gtbone_cyclepart_values, [max_length 1], 'bicubic');
    selectDoF_cyclepart_gt_valuesstretched        = imresize(selectDoF_cyclepart_gt_values, [max_length 1], 'bicubic');

    % store the value
    rmse_us2regbone_cyclepartsnew(:,i)        = rmse_us2regbone_cyclepart_valuesstretched;
    rmse_regbone2gtbone_cyclepartsnew(:,i)    = rmse_regbone2gtbone_cyclepart_valuesstretched;
    selectDoF_kneeJoint_gt_cyclepartsnew(:,i) = selectDoF_cyclepart_gt_valuesstretched;
end



% calculate mean and std for registration error
rmse_us2regbone_cyclemean = mean(rmse_us2regbone_cyclepartsnew, 2);
rmse_us2regbone_cyclestd  = std(rmse_us2regbone_cyclepartsnew, [], 2);

% calculate mean and std for registration error
rmse_regbone2gtbone_cyclemean = mean(rmse_regbone2gtbone_cyclepartsnew, 2);
rmse_regbone2gtbone_cyclestd  = std(rmse_regbone2gtbone_cyclepartsnew, [], 2);

% calculate mean and std for the cycle
selectDoF_kneeJoint_gt_cyclemean = mean(selectDoF_kneeJoint_gt_cyclepartsnew, 2);
selectDoF_kneeJoint_gt_cyclestd  = std(selectDoF_kneeJoint_gt_cyclepartsnew, [], 2);


% display
yyaxis(ax1(1),'left')
% display_shadedError(ax1(1), 1:max_length, selectDoF_kneeJoint_gt_cyclemean, selectDoF_kneeJoint_gt_cyclestd, 'colors', [0 0 1]);
plot(ax1(1), selectDoF_kneeJoint_gt_cyclemean, '-', 'Color', 'b', 'LineWidth', 2);

yyaxis(ax1(1),'right')
bar(ax1(1), rmse_us2regbone_cyclemean, 'BarWidth', 1, 'EdgeColor','none', 'FaceAlpha', 0.2);
bar(ax1(1), rmse_regbone2gtbone_cyclemean, 'BarWidth', 1, 'EdgeColor','none', 'FaceAlpha', 0.2);
ylim(ax1(1), [0, 5]);