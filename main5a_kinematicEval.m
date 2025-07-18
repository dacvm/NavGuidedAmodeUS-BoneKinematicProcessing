%% HELOOOOOO
% - Before running this script, you must already generated .mat file called
%   "all_kneeJoint6DOFs_sX_mXX.mat". This .mat file contains transformations 
%   of the bones (femur and tibia) relative to the ref (global) both GT 
%   and est.
% - The purpose of this script is to quantify and evaluate the kinematic 
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
% dir_Tdata     = 'Tdata_s4_m04_20250717-154844';

% [EDIT] for saving the resulting mat file
is_saveMat = false;

% [EDIT]
is_split  = false;

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

% 1) Load all the data related to knee joint 6 dof estimation.
% -> This mat file is generated from main4_kinematicEstimation.m
tmp_str = split(dir_trial, '_');
sess_str = tmp_str{3};
meas_str = tmp_str{4};
mat_filename = sprintf('all_kneeJoint6DOFs_s%s_m%s.mat', sess_str(end), meas_str);
mat_fullpath = fullfile(path_outputdepth, dir_Tdata, mat_filename);
load(mat_fullpath);

% 2) Load the cycle timestamp data. 
% -> This data indicates where the cycle motion starts, what was been 
%    done in the experiment.
% -> this is obtained from extra_detectFErotCycle.m
csv_filename = sprintf('cycle_timestamp_s%s_m%s.csv', sess_str(end), meas_str);
csv_fullpath = fullfile(path_outputs, csv_filename);
cycle_timestamp = readmatrix(csv_fullpath);
% cycle_timestamp = [158,	353,	543,	749,	957,	1158,	1362,	1565,	1755,	1926,	2086,	2243];


% 3) Alternative color scheme for some plots
c =  [116, 185, 255;
      255, 118, 117;
      253, 203, 110;
       85, 239, 196;
      253, 121, 168;
      162, 155, 254;
      178, 190, 195;
      129, 236, 236] / 255;  

%% INITIALIZE FIGURE OBJECTS AND EVERYTHING RELATED TO PLOTS

% Title for each plot, can be used in any figure
ax_title = {'Medial-Lateral', 'Anterior-Posterior', 'Distraction-Compression', 'Flexion-Extension', 'Abduction-Adduction', 'External-Internal'};

% First figure is to show the joint kinematic with all of the cycle parts
fig1 = figure('Name', 'Joint Kinematic: All Cycle Parts');
t1 = tiledlayout(fig1, 2, 3, ...
     'TileSpacing', 'compact', ...   % tighten spacing if you like
     'Padding',     'compact');      % remove outer margins

% Pre-allocate an array of axes handles
ax1 = gobjects(6,1);

% Populate each tile and store its axes handle
for idx_dof = 1:6
    ax1(idx_dof) = nexttile(t1, idx_dof);
    hold(ax1(idx_dof), 'on');
    grid(ax1(idx_dof), 'on');
    axis(ax1(idx_dof), 'tight');
    title(ax1(idx_dof), ax_title{idx_dof});
    xlabel(ax1(idx_dof), 'Timestamp');
    if(idx_dof<=3)
        ylabel(ax1(idx_dof), 'mm');
    else
        ylabel(ax1(idx_dof), 'deg');
    end
    ax1(idx_dof).FontSize = 12;
end


% Second figure is to show the joint kinematic with the mean and the std of
% all cycle parts
fig2 = figure('Name', 'Joint Kinematic: Mean and STD');
t2 = tiledlayout(fig2, 2, 3, ...
     'TileSpacing', 'compact', ...   % tighten spacing if you like
     'Padding',     'compact');      % remove outer margins

% Pre-allocate an array of axes handles
ax2 = gobjects(6,1);

% Populate each tile and store its axes handle
for idx_dof = 1:6
    ax2(idx_dof) = nexttile(t2, idx_dof);
    hold(ax2(idx_dof), 'on');
    grid(ax2(idx_dof), 'on');
    axis(ax2(idx_dof), 'tight');
    title(ax2(idx_dof), ax_title{idx_dof});
    xlabel(ax2(idx_dof), 'Timestamp');
    ylim(ax2(idx_dof), [-10, 10]);
    if(idx_dof<=3)
        ylabel(ax2(idx_dof), 'mm');
    else
        ylabel(ax2(idx_dof), 'deg');
    end
    ax2(idx_dof).FontSize = 12;
end



%% MAIN PROGRAM

if(~is_split)

% Grab the valid index
timestamp_idcs_valid = all_kneeJoint6DOFs_table.Timestamp_idx;
timestamp_ms_valid = all_kneeJoint6DOFs_table.Timestamp_ms;

% Grab the knee joint 6dof and convert it into matrix
kneeJoint6DOFs_est = cell2mat(all_kneeJoint6DOFs_table.kneeJoint6DOFs_est);
kneeJoint6DOFs_gt  = cell2mat(all_kneeJoint6DOFs_table.kneeJoint6DOFs_gt);

cycle_select     = [2,5];
cycle_idx_select = [cycle_timestamp(cycle_select(1)), cycle_timestamp(cycle_select(2))];

for idx_dof=1:6

    % get the current dof
    currentDoF_kneeJoint_est = kneeJoint6DOFs_est(cycle_idx_select(1):cycle_idx_select(2), idx_dof);
    currentDoF_kneeJoint_gt = kneeJoint6DOFs_gt(cycle_idx_select(1):cycle_idx_select(2), idx_dof);
    % currentDoF_kneeJoint_est = kneeJoint6DOFs_est(1:500, idx_dof);
    % currentDoF_kneeJoint_gt = kneeJoint6DOFs_gt(1:500, idx_dof);

    % calculate the difference
    currentDoF_kneeJoint_estgt_diff   = currentDoF_kneeJoint_est - currentDoF_kneeJoint_gt;

    % plot the knee joint for the current dof
    plot(ax1(idx_dof), currentDoF_kneeJoint_est, '-', 'Color', 'r', 'LineWidth', 2);
    plot(ax1(idx_dof), currentDoF_kneeJoint_gt, '-', 'Color', 'b', 'LineWidth', 2);
    legend(ax1(idx_dof), {'Estimation', 'Ground truth'});

    % compute the quantitative result
    result_corr  = corr(currentDoF_kneeJoint_est, currentDoF_kneeJoint_gt);
    result_rmse  = sqrt( mean( (currentDoF_kneeJoint_estgt_diff).^2 ) );
    result_std   = std(currentDoF_kneeJoint_estgt_diff);
    [~, ~, ~, result_q2, ~, result_uw, ~] = computeBoxplotStats(abs(currentDoF_kneeJoint_estgt_diff));

    str = sprintf('corr   = %.2f\nrmse = %.2f ± %.2f\nmad  = %.2f%s%.2f', result_corr, result_rmse, result_std, result_q2, char(9652), result_uw);
    text( ax1(idx_dof), ...
          0.05, 0.95, ... 
          str, ...
          'Units',            'normalized', ...
          'VerticalAlignment','top', ...
          'FontSize',         10, ...
          'BackgroundColor',  'w', ...
          'EdgeColor',        'k', ...
          'Margin',           5 );

    str = sprintf('rmse = %.2f ± %.2f\nmad  = %.2f%s%.2f', result_rmse, result_std, result_q2, char(9652), result_uw);
    text( ax2(idx_dof), ...
          0.05, 0.95, ... 
          str, ...
          'Units',            'normalized', ...
          'VerticalAlignment','top', ...
          'FontSize',         10, ...
          'BackgroundColor',  'w', ...
          'EdgeColor',        'k', ...
          'Margin',           5 );

end

% we don't need the fig2 if the user choose the "non-split" plot
close(fig2);


else
%% Figure 1, kinematic plot

% Grab the valid index
timestamp_idcs_valid = all_kneeJoint6DOFs_table.Timestamp_idx;

% Grab the knee joint 6dof and convert it into matrix
kneeJoint6DOFs_est = cell2mat(all_kneeJoint6DOFs_table.kneeJoint6DOFs_est);
kneeJoint6DOFs_gt  = cell2mat(all_kneeJoint6DOFs_table.kneeJoint6DOFs_gt);

% Select cycle
cycle_select      = [2,5];
cycle_idcs_select = cycle_timestamp(cycle_select(1):cycle_select(2));

for idx_dof=1:6

    % get the current dof
    currentDoF_kneeJoint_est = kneeJoint6DOFs_est(:,idx_dof);
    currentDoF_kneeJoint_gt = kneeJoint6DOFs_gt(:,idx_dof);

    % before we do further processing, i want to chop the knee joint data
    % into different cyclic motion parts. here, i initialize a cell to 
    % store all the parts (parts can have different length)
    currentDoF_kneeJoint_est_cycleparts = {};
    currentDoF_kneeJoint_gt_cycleparts  = {};

    % get the cycle timestamp
    for i=1:(length(cycle_idcs_select)-1)

        % get the start and end index
        idx_start = cycle_idcs_select(i);
        idx_end = cycle_idcs_select(i+1);
        
        % before we go further let's check if the end index is bigger than
        % the valid timestamp we have
        if(idx_end>timestamp_idcs_valid(end))
            % if yes, let's not use this specific cycle part, stop it.
            break;
        end

        % get the current indices of the cycle part
        currentDoF_cyclepart_idcs  = find( (timestamp_idcs_valid>=idx_start) & (timestamp_idcs_valid<idx_end) );
        % get the current value of the cycle part
        currentDoF_cyclepart_est_values = currentDoF_kneeJoint_est(currentDoF_cyclepart_idcs);
        currentDoF_cyclepart_gt_values  = currentDoF_kneeJoint_gt(currentDoF_cyclepart_idcs);

        % store the current values of this cycle part
        currentDoF_kneeJoint_est_cycleparts{i} = currentDoF_cyclepart_est_values;
        currentDoF_kneeJoint_gt_cycleparts{i}  = currentDoF_cyclepart_gt_values;
    end

    % get how many parts we have (both gt and est must be the same)
    n_cycleparts = length(currentDoF_kneeJoint_est_cycleparts);
    % get the max length of the vector inside the cell
    max_length = max( cellfun(@numel, currentDoF_kneeJoint_est_cycleparts) );

    % allocate memory to store the new stretched version of the part
    currentDoF_kneeJoint_est_cyclepartsnew = zeros(max_length, n_cycleparts);
    currentDoF_kneeJoint_gt_cyclepartsnew = zeros(max_length, n_cycleparts);

    % we will loop for each cycle parts and we will stretch the shorter 
    % parts to be the same length as the longest part.
    for i=1:n_cycleparts

        % get the current cycle part
        currentDoF_cyclepart_est_values = currentDoF_kneeJoint_est_cycleparts{i};
        currentDoF_cyclepart_gt_values = currentDoF_kneeJoint_gt_cycleparts{i};
        
        % stretch the value, here i am using imresize from image processing
        % tool box, because why not? they already have built-in
        % interpolation function inside
        currentDoF_cyclepart_est_valuesstretched = imresize(currentDoF_cyclepart_est_values, [max_length 1], 'bicubic');
        currentDoF_cyclepart_gt_valuesstretched  = imresize(currentDoF_cyclepart_gt_values, [max_length 1], 'bicubic');

        % store the value
        currentDoF_kneeJoint_est_cyclepartsnew(:,i) = currentDoF_cyclepart_est_valuesstretched;
        currentDoF_kneeJoint_gt_cyclepartsnew(:,i)  = currentDoF_cyclepart_gt_valuesstretched;

        % plot the stretched version
        plot(ax1(idx_dof), currentDoF_cyclepart_est_valuesstretched, '-', 'Color', '#FFCDD2', 'LineWidth', 0.5);
        plot(ax1(idx_dof), currentDoF_cyclepart_gt_valuesstretched, '-', 'Color', '#BBDEFB', 'LineWidth', 0.5);
    end

    % calculate mean and std for the estimation...
    currentDoF_kneeJoint_est_cyclemean = mean(currentDoF_kneeJoint_est_cyclepartsnew, 2);
    currentDoF_kneeJoint_est_cyclestd  = std(currentDoF_kneeJoint_est_cyclepartsnew, [], 2);
    % ...and for the gt
    currentDoF_kneeJoint_gt_cyclemean = mean(currentDoF_kneeJoint_gt_cyclepartsnew, 2);
    currentDoF_kneeJoint_gt_cyclestd  = std(currentDoF_kneeJoint_gt_cyclepartsnew, [], 2);
    % then calculate the difference between the mean
    currentDoF_kneeJoint_estgt_diff   = currentDoF_kneeJoint_est_cyclemean - currentDoF_kneeJoint_gt_cyclemean;

    % plot the mean and std
    % display_shadedError(ax1(idx_dof), 1:max_length, current_kneeJointDOF_est_cyclemean, current_kneeJointDOF_est_cyclestd, 'colors', [1 0 0]);
    % display_shadedError(ax1(idx_dof), 1:max_length, current_kneeJointDOF_gt_cyclemean, current_kneeJointDOF_gt_cyclestd, 'colors', [0 0 1]);
    plot(ax1(idx_dof), currentDoF_kneeJoint_est_cyclemean, '-', 'Color', 'r', 'LineWidth', 2);
    plot(ax1(idx_dof), currentDoF_kneeJoint_gt_cyclemean, '-', 'Color', 'b', 'LineWidth', 2);
    % set the legend, but i only want the last two (stored with newest first)
    mylines = findobj(ax1(idx_dof), 'Type','Line');
    legend(mylines(1:2), {'Ground truth','Estimation'}, 'FontSize', 12);

    % plot the error, relative to the ground truth
    plot(ax2(idx_dof), currentDoF_kneeJoint_estgt_diff, '-', 'Color', 'r', 'LineWidth', 2);
    plot(ax2(idx_dof), zeros(max_length, 1), '-', 'Color', 'b', 'LineWidth', 2);
    % set the legend, but i only want the last two (stored with newest first)
    mylines = findobj(ax2(idx_dof), 'Type','Line');
    legend(mylines(1:2), {'Ground truth','Estimation'});

    % compute the quantitative result
    result_corr  = corr(currentDoF_kneeJoint_est_cyclemean, currentDoF_kneeJoint_gt_cyclemean);
    result_rmse  = sqrt( mean( (currentDoF_kneeJoint_estgt_diff).^2 ) );
    result_std   = std(currentDoF_kneeJoint_estgt_diff);
    [~, ~, ~, result_q2, ~, result_uw, ~] = computeBoxplotStats(abs(currentDoF_kneeJoint_estgt_diff));

    str = sprintf('corr   = %.2f\nrmse = %.2f ± %.2f\nmad  = %.2f%s%.2f', result_corr, result_rmse, result_std, result_q2, char(9652), result_uw);
    text( ax1(idx_dof), ...
          0.05, 0.95, ... 
          str, ...
          'Units',            'normalized', ...
          'VerticalAlignment','top', ...
          'FontSize',         10, ...
          'BackgroundColor',  'w', ...
          'EdgeColor',        'k', ...
          'Margin',           5 );

    str = sprintf('rmse = %.2f ± %.2f\nmad  = %.2f%s%.2f', result_rmse, result_std, result_q2, char(9652), result_uw);
    text( ax2(idx_dof), ...
          0.05, 0.95, ... 
          str, ...
          'Units',            'normalized', ...
          'VerticalAlignment','top', ...
          'FontSize',         10, ...
          'BackgroundColor',  'w', ...
          'EdgeColor',        'k', ...
          'Margin',           5 );


    
end


% end if(~is_split)
end


% activate the figure to be on top
figure(fig1);