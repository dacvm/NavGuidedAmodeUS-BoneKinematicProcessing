%% HELOOOOOO
% - Before running this script, you must already generated .mat files:
%   1) all_kneeJoint6DOFs_sX_mXX.mat
%      - Contains transformations of the bones (femur and tibia) relative to 
%        the ref (global), and the 6Dof of knee joint kinematic. 
%      - Generated from main4_kinematicEstimation3_TFJointRUMC.m
%   2) cycle_timestamp_sX_mXX.csv
%      - Contains the timestamps of the cycle. 
%       - Generated from extra_detectFErotCycle.m
% - This script will produce 6 plots, each plot for each individual degree
%   of freedom of knee joint kinematic estimation against its ground truth
% - There are two options to show the plot: 
%   1) all continuous cycle shown at once, 
%   2) Cycles are chopped and each chop will be shown overalayed to each
%   other together with their mean.
% - If you want to generate the error relative to ground truth, check
%   main5a_kinematicEvalAll2.m

clc; clear; close all;

% [EDIT] directory to the project
path_root    = 'D:\Documents\BELANDA\PhD Thesis\Code\MATLAB\amode_navigation_experiment\experiment_b';

% [EDIT] Change the data you are using accordingly. 
% -----> dir_depthdata is created by main1_processDepthData.m
% -----> dir_Tdata is created by main3_registrationWithTime.m
dirs_Tdata = { fullfile('depthdata_s4_m04_20250708-172830', 'Tdata_s4_m04_20250808-212946'), ...   % with-nav
               fullfile('depthdata_s3_m02_20250722-114731', 'Tdata_s3_m02_20250801-203337'), ...   % no-nav, manual
               fullfile('depthdata_s3_m02_20250722-174503', 'Tdata_s3_m02_20250807-215603'), ...   % no-nav, auto, 2x noise
               fullfile('depthdata_s3_m02_20250722-174812', 'Tdata_s3_m02_20250810-215516')};      % no-nav, auto, 3x noise

% [EDIT] select which data you want to show
idx_dir_Tdata = 1;

% [EDIT] select which cycle you want to show
cycle_select  = [3,12];

% [EDIT]
% do you want to split the data into cycles and average them?
is_split  = false;
% if true, do you want to show the std as a shade or washout-lines?
is_shaded = true;

% [EDIT] select the color scheme
color_scheme = 2;

% [EDIT] for saving the resulting mat file
is_saveMat = false;


%% INITIALIZE PATHS AND LOADING SOME CONFIGURATION

% declare some of the important paths
path_function     = fullfile(path_root, 'functions');
path_outputs      = fullfile(path_root, 'outputs');
path_data         = fullfile(path_root, 'data');

% Generate path to function directory
addpath(genpath(path_function));


%% INITIALIZE THE DATA TO SHOW

% 1) Load all the data related to knee joint 6 dof estimation.
% -> This mat file is generated from main4_kinematicEstimation.m
tmp_str  = strsplit(dirs_Tdata{idx_dir_Tdata}, {'_', '\'});
sess_str = tmp_str{2};
meas_str = tmp_str{3};
mat_filename = sprintf('all_kneeJoint6DOFs_%s_%s.mat', sess_str, meas_str);
mat_fullpath = fullfile(path_outputs, 'output_allest', dirs_Tdata{idx_dir_Tdata}, mat_filename);
load(mat_fullpath);

% 2) Load the cycle timestamp data. 
% -> This data indicates where the cycle motion starts, what was been 
%    done in the experiment.
% -> this is obtained from extra_detectFErotCycle.m
csv_filename = sprintf('cycle_timestamp_%s_%s.csv', sess_str, meas_str);
csv_fullpath = fullfile(path_outputs, csv_filename);
cycle_timestamp = readmatrix(csv_fullpath);


%% INITIALIZE FIGURE OBJECTS AND EVERYTHING RELATED TO PLOTS

% Title for each plot, can be used in any figure. The title will be depend
% on how do we calculate the knee joint kinematic. The kneeJoint_method
% value determine with which method we use, the title will follow
% accordingly
if(kneeJoint_method==1)
    ax_title = { 'Medial - Lateral', ...
                 'Anterior - Posterior', ...
                 'Distraction - Compression', ...
                 'Flexion - Extension', ...
                 'Abduction - Adduction', ...
                 'External-Internal'};
elseif (kneeJoint_method==2)
    ax_title = { 'Anterior(+) - Posterior(-)', ...
                 'Proximal(+) - Distal(-)', ...
                 'Medial(+) - Lateral(-)', ...
                 'Flexion(+) - Extension(-)', ...
                 'Valgus(+) - Varus(-)', ...
                 'Exorotaion(+) - Endorotation(-)'};
elseif (kneeJoint_method==3)
    ax_title = {'Anterior(+) - Posterior(-)', ... 
                'Proximal(+) - Distal(-)', ...
                'Lateral(+) - Medial(-)', ...
                'Flexion(+) - Extension(-)', ...
                'External(+) - Internal(-)', ...
                'Abduction(+) - Adduction(-)'};
end

% color scheme (light)
if(color_scheme==1)
c =  [116, 185, 255;
      255, 118, 117;
      253, 203, 110;
       85, 239, 196;
      253, 121, 168;
      162, 155, 254;
      178, 190, 195;
      129, 236, 236] / 255;
% color scheme (medium)
elseif(color_scheme==2)
c = [ 52, 152, 219;
     231,  76,  60;
     243, 156,  18;
      46, 204, 113;
     155,  89, 182;
      52,  73,  94] / 255;
% color scheme (dark)
else
c = [ 41, 128, 185
     192,  57,  43
     230, 126,  34
      41, 128, 185
     142,  68, 173
      44,  62, 80] /255;
end
    
% set the ylim so that the axes become consistent
ax_ylim = repmat([-20, 20], 6, 1);

% First figure is to show the joint kinematic with all of the cycle parts
fig1 = figure('Name', 'Joint Kinematic: All Cycle Parts', 'Position', [50 300 1700 700]);
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
    xlabel(ax1(idx_dof), 'Time (s)');
    if(idx_dof<=3)
        ylabel(ax1(idx_dof), 'mm');
    else
        ylabel(ax1(idx_dof), 'deg');
    end
    ax1(idx_dof).FontSize = 12;
    % ylim( ax1(idx_dof), ax_ylim(idx_dof,:) );
end


%% MAIN PROGRAM
% I tried to make the numbering (steps of program) consistent with
% main5a_kinematicEvalAll2.m. Hopefuly by doing this it is easier to
% read this script and that script.
% But in all honesty, i am focus more on this non-splitting part, as it is
% more logical to evaluate the raw data instead of average (splitting).
% So forgive me if the splitting part is a bit confusing and inconsistent
% with the other scripts.


% 1) Select cycle (In this part cycle_timestamp already loaded)
cycle_idxoriginal_select = [cycle_timestamp(cycle_select(1)), cycle_timestamp(cycle_select(2))];

% 2) Delete problematic rows (In this part all_kneeJoint6DOFs_table already loaded.)
idcs_problematicRows = find(all_kneeJoint6DOFs_table.is_invalid);
all_kneeJoint6DOFs_table(idcs_problematicRows, :) = [];

% To store the metric evaluation values.
% rows will be organized as: tx, ty, tz, rx, ry, rz 
% (be aware, i choose xyz here, easier for reader to understand, 
% instead of the default format, zyx, from matlab)
% cols will be organized as: rmse, std, mad, uw, corr
metric_mat = zeros(6, 5);


% If the user select splitting
if(~is_split)

% 3) Grab the valid timestamp
timestamp_idcs_valid = all_kneeJoint6DOFs_table.Timestamp_idx;
timestamp_ms_valid   = all_kneeJoint6DOFs_table.Timestamp_ms;


% 4) The indices that are shown in cycle_timestamps are the indices of the 
%    table of the original measurement (purely from experiment without
%    any processing at all). 
% -> We make a lot of filtering process already, for example, removing some 
%    invalid data. So, the all_kneeJoint6DOFs_table.Timestamp_idx might
%    have holes within it.
% -> By performing the process below, we will get to get the actual timestamp.
tmp1 = find(timestamp_idcs_valid==cycle_idxoriginal_select(1));
tmp2 = find(timestamp_idcs_valid==cycle_idxoriginal_select(2));
cycle_idxvalid_select = [tmp1, tmp2];

% 5) Grab the knee joint 6dof and convert it into matrix
kneeJoint6DOFs_est = cell2mat(all_kneeJoint6DOFs_table.kneeJoint6DOFs_est);
kneeJoint6DOFs_gt  = cell2mat(all_kneeJoint6DOFs_table.kneeJoint6DOFs_gt);

% 6) Get the selected timestamp_ms
timestamp_ms_valid_select  = timestamp_ms_valid(cycle_idxvalid_select(1):cycle_idxvalid_select(2));
timestamp_ms_valid_select0 = timestamp_ms_valid_select - timestamp_ms_valid_select(1);

% loop for all dofs, let's display them one by one
for idx_dof=1:6

    % 1) Get the current dof
    currentDoF_kneeJoint_est = kneeJoint6DOFs_est(cycle_idxvalid_select(1):cycle_idxvalid_select(2), idx_dof);
    currentDoF_kneeJoint_est = smoothdata(currentDoF_kneeJoint_est, 'rlowess', 20);

    % 2) Get the current gt dof
    currentDoF_kneeJoint_gt = kneeJoint6DOFs_gt(cycle_idxvalid_select(1):cycle_idxvalid_select(2), idx_dof);

    % 3) Calculate the difference
    currentDoF_kneeJoint_estgt_diff   = currentDoF_kneeJoint_est - currentDoF_kneeJoint_gt;

    % 4) Plot the knee joint for the current dof
    plot(ax1(idx_dof), timestamp_ms_valid_select0*1e-3, currentDoF_kneeJoint_est, '-', 'Color', 'r', 'LineWidth', 2);
    plot(ax1(idx_dof), timestamp_ms_valid_select0*1e-3, currentDoF_kneeJoint_gt, '-', 'Color', 'b', 'LineWidth', 2);
    legend(ax1(idx_dof), {'Estimation', 'Ground truth'});

    % compute the quantitative result
    result_corr  = corr(currentDoF_kneeJoint_est, currentDoF_kneeJoint_gt);
    result_rmse  = sqrt( mean( (currentDoF_kneeJoint_estgt_diff).^2 ) );
    result_std   = std(currentDoF_kneeJoint_estgt_diff);
    [~, ~, ~, result_q2, ~, result_uw, ~] = computeBoxplotStats(abs(currentDoF_kneeJoint_estgt_diff));

    str = sprintf('corr   = %.2f\nrmse = %.2f %s %.2f\nmad  = %.2f%s%.2f', result_corr, result_rmse, char(177), result_std, result_q2, char(9652), result_uw);
    text( ax1(idx_dof), ...
          0.05, 0.95, ... 
          str, ...
          'Units',            'normalized', ...
          'VerticalAlignment','top', ...
          'FontSize',         10, ...
          'BackgroundColor',  'w', ...
          'EdgeColor',        'k', ...
          'Margin',           5 );

    % store the value
    metric_mat(idx_dof, :) = [result_rmse, result_std, result_q2, result_uw, result_corr]';

end


else
%% Figure 1, kinematic plot

% Grab the valid index
timestamp_idcs_valid = all_kneeJoint6DOFs_table.Timestamp_idx;
timestamp_ms_valid   = all_kneeJoint6DOFs_table.Timestamp_ms;

% Grab the knee joint 6dof and convert it into matrix
kneeJoint6DOFs_est = cell2mat(all_kneeJoint6DOFs_table.kneeJoint6DOFs_est);
kneeJoint6DOFs_gt  = cell2mat(all_kneeJoint6DOFs_table.kneeJoint6DOFs_gt);

% Select cycle
cycle_idxoriginal_select = cycle_timestamp(cycle_select(1):cycle_select(2));

for idx_dof=1:6

    % get the current dof
    currentDoF_kneeJoint_est = kneeJoint6DOFs_est(:,idx_dof);
    currentDoF_kneeJoint_gt  = kneeJoint6DOFs_gt(:,idx_dof);

    % before we do further processing, i want to chop the knee joint data
    % into different cyclic motion parts. here, i initialize a cell to 
    % store all the parts (parts can have different length)
    currentDoF_kneeJoint_est_cycleparts = {};
    currentDoF_kneeJoint_gt_cycleparts  = {};
    timestamp_ms_cycleparts             = {};

    % get the cycle timestamp
    for i=1:(length(cycle_idxoriginal_select)-1)

        % Get the start and end index
        % The indices that are shown in cycle_timestamp is the indices of 
        % the table original measurement. We make a lot of filtering 
        % process already (removing some invalid data). 
        % So we need to get the actual timestamp.
        idx_start = find(timestamp_idcs_valid==cycle_idxoriginal_select(i));
        idx_end = find(timestamp_idcs_valid==cycle_idxoriginal_select(i+1));

        % % get the start and end index
        % idx_start = cycle_idxvalid_select(i);
        % idx_end = cycle_idxvalid_select(i+1);
        
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
        timestamp_cyclepart_ms_values   = timestamp_ms_valid(currentDoF_cyclepart_idcs);

        % store the current values of this cycle part
        currentDoF_kneeJoint_est_cycleparts{i} = currentDoF_cyclepart_est_values;
        currentDoF_kneeJoint_gt_cycleparts{i}  = currentDoF_cyclepart_gt_values;
        timestamp_ms_cycleparts{i}             = timestamp_cyclepart_ms_values;
    end

    % get how many parts we have (both gt and est must be the same)
    n_cycleparts = length(currentDoF_kneeJoint_est_cycleparts);
    % get the max length of the vector inside the cell
    [max_length_val, max_length_idx] = max( cellfun(@numel, currentDoF_kneeJoint_est_cycleparts) );

    % get the longest timestamp, this vector will be used for x-axis for plot
    timestamp_ms_maxlength  = timestamp_ms_cycleparts{max_length_idx};
    timestamp_ms_maxlength0 = timestamp_ms_maxlength - timestamp_ms_maxlength(1);

    % allocate memory to store the new stretched version of the part
    currentDoF_kneeJoint_est_cyclepartsnew = zeros(max_length_val, n_cycleparts);
    currentDoF_kneeJoint_gt_cyclepartsnew = zeros(max_length_val, n_cycleparts);

    % we will loop for each cycle parts and we will stretch the shorter 
    % parts to be the same length as the longest part.
    for i=1:n_cycleparts

        % get the current cycle part
        currentDoF_cyclepart_est_values = currentDoF_kneeJoint_est_cycleparts{i};
        currentDoF_cyclepart_gt_values = currentDoF_kneeJoint_gt_cycleparts{i};
        
        % stretch the value, here i am using imresize from image processing
        % tool box, because why not? they already have built-in
        % interpolation function inside
        currentDoF_cyclepart_est_valuesstretched = imresize(currentDoF_cyclepart_est_values, [max_length_val 1], 'bicubic');
        currentDoF_cyclepart_gt_valuesstretched  = imresize(currentDoF_cyclepart_gt_values, [max_length_val 1], 'bicubic');

        % store the value
        currentDoF_kneeJoint_est_cyclepartsnew(:,i) = currentDoF_cyclepart_est_valuesstretched;
        currentDoF_kneeJoint_gt_cyclepartsnew(:,i)  = currentDoF_cyclepart_gt_valuesstretched;

        % plot the stretched version
        if(~is_shaded)
            plot(ax1(idx_dof), timestamp_ms_maxlength0 * 1e-3, currentDoF_cyclepart_est_valuesstretched, '-', 'Color', '#FFEBEE', 'LineWidth', 0.5);
            plot(ax1(idx_dof), timestamp_ms_maxlength0 * 1e-3, currentDoF_cyclepart_gt_valuesstretched, '-', 'Color', '#E3F2FD', 'LineWidth', 0.5);
        end
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
    if(is_shaded)
        display_shadedError(ax1(idx_dof), timestamp_ms_maxlength0 * 1e-3, currentDoF_kneeJoint_est_cyclemean, currentDoF_kneeJoint_est_cyclestd, 'colors', [1 0 0]);
        display_shadedError(ax1(idx_dof), timestamp_ms_maxlength0 * 1e-3, currentDoF_kneeJoint_gt_cyclemean, currentDoF_kneeJoint_gt_cyclestd, 'colors', [0 0 1]);
    else
        plot(ax1(idx_dof), timestamp_ms_maxlength0 * 1e-3, currentDoF_kneeJoint_est_cyclemean, '-', 'Color', 'r', 'LineWidth', 2);
        plot(ax1(idx_dof), timestamp_ms_maxlength0 * 1e-3, currentDoF_kneeJoint_gt_cyclemean, '-', 'Color', 'b', 'LineWidth', 2);
    end
    % set the legend, but i only want the last two (stored with newest first)
    mylines = findobj(ax1(idx_dof), 'Type','Line');
    legend(mylines(1:2), {'Ground truth','Estimation'}, 'FontSize', 12);

    % compute the quantitative result
    result_corr  = corr(currentDoF_kneeJoint_est_cyclemean, currentDoF_kneeJoint_gt_cyclemean);
    result_rmse  = sqrt( mean( (currentDoF_kneeJoint_estgt_diff).^2 ) );
    result_std   = std(currentDoF_kneeJoint_estgt_diff);
    [~, ~, ~, result_q2, ~, result_uw, ~] = computeBoxplotStats(abs(currentDoF_kneeJoint_estgt_diff));

    str = sprintf('corr   = %.2f\nrmse = %.2f %s %.2f\nmad  = %.2f%s%.2f', result_corr, result_rmse, char(177), result_std, result_q2, char(9652), result_uw);
    text( ax1(idx_dof), ...
          0.05, 0.95, ... 
          str, ...
          'Units',            'normalized', ...
          'VerticalAlignment','top', ...
          'FontSize',         10, ...
          'BackgroundColor',  'w', ...
          'EdgeColor',        'k', ...
          'Margin',           5 );
    
    % store the value
    metric_mat(idx_dof, :) = [result_rmse, result_std, result_q2, result_uw, result_corr]';
end


% end if(~is_split)
end

% Swap row-4 and row-6, from r-zyx to r-xyz
tmp = metric_mat(4, :);
metric_mat(4, :) = metric_mat(6, :);
metric_mat(6, :) = tmp;

% Create a table. You can just open the table and copy paste to excel 
% or latex table generator.
metric_table = array2table( metric_mat,  'VariableNames', ...
                            {'rmse_mean', 'rmse_std', 'mad_q2', 'mad_uw', 'corr'});