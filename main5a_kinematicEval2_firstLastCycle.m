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
cycle_select  = [3, 6; 9, 12];

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



% 1) Delete problematic rows (In this part all_kneeJoint6DOFs_table already loaded.)
idcs_problematicRows = find(all_kneeJoint6DOFs_table.is_invalid);
all_kneeJoint6DOFs_table(idcs_problematicRows, :) = [];

% 2) Grab the valid timestamp
timestamp_idcs_valid = all_kneeJoint6DOFs_table.Timestamp_idx;
timestamp_ms_valid   = all_kneeJoint6DOFs_table.Timestamp_ms;

% 3) Grab the knee joint 6dof and convert it into matrix
kneeJoint6DOFs_est = cell2mat(all_kneeJoint6DOFs_table.kneeJoint6DOFs_est);
kneeJoint6DOFs_gt  = cell2mat(all_kneeJoint6DOFs_table.kneeJoint6DOFs_gt);




% 4) Select cycle (In this part cycle_timestamp already loaded)
cycle_idxoriginal_select = [ cycle_timestamp(cycle_select(1, 1)), ...
                             cycle_timestamp(cycle_select(1, 2)) ];
% 5) The indices that are shown in cycle_timestamps (and also in 
%    cycle_idxoriginal_select) are the indices of the table of the 
%    original measurement (purely from experiment without any 
%    processing at all). 
% -> We make a lot of filtering process already, for example:
%    - In the depth detection part, we deleted some parts that are an 
%      artefact from the extrapolation function (shoot up and down). 
%    - Just before this, we remove some invalid data. So, the 
%      all_kneeJoint6DOFs_table.Timestamp_idx might have holes within it.
% -> By performing the process below, we will get the table index based
%    from the value shown in the cycle_timestamps
tmp1 = find(timestamp_idcs_valid==cycle_idxoriginal_select(1));
tmp2 = find(timestamp_idcs_valid==cycle_idxoriginal_select(2));
cycle_idxvalid_select_part1 = [tmp1, tmp2];

% Repeat step (4)
cycle_idxoriginal_select = [ cycle_timestamp(cycle_select(2, 1)), ...
                             cycle_timestamp(cycle_select(2, 2)) ];
% Repeat step (5)
tmp1 = find(timestamp_idcs_valid==cycle_idxoriginal_select(1));
tmp2 = find(timestamp_idcs_valid==cycle_idxoriginal_select(2));
cycle_idxvalid_select_part2 = [tmp1, tmp2];


% To store the metric evaluation values.
% rows will be organized as: tx, ty, tz, rx, ry, rz 
% (be aware, i choose xyz here, easier for reader to understand, 
% instead of the default format, zyx, from matlab)
% cols will be organized as: rmse, std, mad, uw, corr
metric_mat_part1 = zeros(6, 5);
metric_mat_part2 = zeros(6, 5);
pvalues_rmse     = zeros(6, 2);
pvalues_mad      = zeros(6, 2);

% loop for all dofs, let's display them one by one
for idx_dof=1:6

        % 1) Get the part 1 -----------------------------------------------

        % 1.1) Get the current dof
        currentDoF_kneeJoint_est_part1 = kneeJoint6DOFs_est(cycle_idxvalid_select_part1(1):cycle_idxvalid_select_part1(2), idx_dof);
        currentDoF_kneeJoint_est_part1 = smoothdata(currentDoF_kneeJoint_est_part1, 'rlowess', 20);

        % 1.2) Get the current gt
        currentDoF_kneeJoint_gt_part1  = kneeJoint6DOFs_gt(cycle_idxvalid_select_part1(1):cycle_idxvalid_select_part1(2), idx_dof);
    
        % 1.3) Get the selected timestamp_ms
        timestamp_ms_valid_select_part1  = timestamp_ms_valid(cycle_idxvalid_select_part1(1):cycle_idxvalid_select_part1(2));
        timestamp_ms_valid_select0_part1 = timestamp_ms_valid_select_part1 - cycle_idxvalid_select_part1(1);
        
        % 1.4) Compute quantitative result
        currentDoF_kneeJoint_estgt_diff1   = currentDoF_kneeJoint_est_part1 - currentDoF_kneeJoint_gt_part1;
        
        % 1.5) Compute the quantitative result
        result_corr_1  = corr(currentDoF_kneeJoint_est_part1, currentDoF_kneeJoint_gt_part1);
        result_rmse_1  = sqrt( mean( (currentDoF_kneeJoint_estgt_diff1).^2 ) );
        result_std_1   = std(currentDoF_kneeJoint_estgt_diff1);
        [~, ~, ~, result_q2_1, ~, result_uw_1, ~] = computeBoxplotStats(abs(currentDoF_kneeJoint_estgt_diff1));

        % 1.6) Store the value
        metric_mat_part1(idx_dof, :) = [result_rmse_1, result_std_1, result_q2_1, result_uw_1, result_corr_1]';


        % 2) Get the part 2 -----------------------------------------------
        % 2.1) Get the current dof
        currentDoF_kneeJoint_est_part2 = kneeJoint6DOFs_est(cycle_idxvalid_select_part2(1):cycle_idxvalid_select_part2(2), idx_dof);
        currentDoF_kneeJoint_est_part2 = smoothdata(currentDoF_kneeJoint_est_part2, 'rlowess', 20);

        % 2.2) Get the current gt
        currentDoF_kneeJoint_gt_part2  = kneeJoint6DOFs_gt(cycle_idxvalid_select_part2(1):cycle_idxvalid_select_part2(2), idx_dof);
    
        % 2.3) Get the selected timestamp_ms
        timestamp_ms_valid_select_part2  = timestamp_ms_valid(cycle_idxvalid_select_part2(1):cycle_idxvalid_select_part2(2));
        timestamp_ms_valid_select0_part2 = timestamp_ms_valid_select_part2 - cycle_idxvalid_select_part1(1);

        % 2.4) Calculate the difference
        currentDoF_kneeJoint_estgt_diff2   = currentDoF_kneeJoint_est_part2 - currentDoF_kneeJoint_gt_part2;

        % 2.5) Compute the quantitative result
        result_corr_2  = corr(currentDoF_kneeJoint_est_part2, currentDoF_kneeJoint_gt_part2);
        result_rmse_2  = sqrt( mean( (currentDoF_kneeJoint_estgt_diff2).^2 ) );
        result_std_2   = std(currentDoF_kneeJoint_estgt_diff2);
        [~, ~, ~, result_q2_2, ~, result_uw_2, ~] = computeBoxplotStats(abs(currentDoF_kneeJoint_estgt_diff2));

        % 2.6) Store the value
        metric_mat_part2(idx_dof, :) = [result_rmse_2, result_std_2, result_q2_2, result_uw_2, result_corr_2]';


        % 3) Display ------------------------------------------------------

        % 3.1a) Display part 1
        plotDiscontinuousXAxis( ax1(idx_dof), ...
                                   timestamp_ms_valid_select0_part1*1e-3, currentDoF_kneeJoint_est_part1, ...
                                   timestamp_ms_valid_select0_part2*1e-3, currentDoF_kneeJoint_est_part2, ...
                                   'GapFraction', 0.03, ...
                                   'NumTicks', 10, ... 
                                   'LineWidth', 2, ...
                                   'Color', [1 0 0], ...
                                   'LineStyle', '-', ...
                                   'ShadeAlpha', 0.05);

        % 3.1b) Display the statistics
        str = sprintf('corr   = %.2f\nrmse = %.2f %s %.2f\nmad  = %.2f%s%.2f', result_corr_1, result_rmse_1, char(177), result_std_1, result_q2_1, char(9652), result_uw_1);
        text( ax1(idx_dof), ...
              0.05, 0.95, ... 
              str, ...
              'Units',            'normalized', ...
              'VerticalAlignment','top', ...
              'FontSize',         10, ...
              'BackgroundColor',  'w', ...
              'EdgeColor',        'k', ...
              'Margin',           5 );

        % 3.2b) Display part 2
        plotDiscontinuousXAxis( ax1(idx_dof), ...
                                   timestamp_ms_valid_select0_part1*1e-3, currentDoF_kneeJoint_gt_part1, ...
                                   timestamp_ms_valid_select0_part2*1e-3, currentDoF_kneeJoint_gt_part2, ...
                                   'GapFraction', 0.03, ...
                                   'NumTicks', 10, ...
                                   'LineWidth', 2, ...
                                   'Color', [0 0 1], ...
                                   'LineStyle', '-', ...
                                   'ShadeAlpha', 0.05);

        % 3.2b) Display the statistics
        str = sprintf('corr   = %.2f\nrmse = %.2f %s %.2f\nmad  = %.2f%s%.2f', result_corr_2, result_rmse_2, char(177), result_std_2, result_q2_2, char(9652), result_uw_2);
        text( ax1(idx_dof), ...
              0.95, 0.95, ... 
              str, ...
              'Units',               'normalized', ...
              'VerticalAlignment',   'top', ...
              'HorizontalAlignment', 'right', ...
              'FontSize',            10, ...
              'BackgroundColor',     'w', ...
              'EdgeColor',           'k', ...
              'Margin',              5 );
        

        % 3.2) Adjust the ylim again
        ylim(ax1(idx_dof), [min([currentDoF_kneeJoint_est_part1; ...
                                 currentDoF_kneeJoint_est_part2; ...
                                 currentDoF_kneeJoint_gt_part1; ...
                                 currentDoF_kneeJoint_gt_part2]), ...
                            max([currentDoF_kneeJoint_est_part1; ...
                                 currentDoF_kneeJoint_est_part2; ...
                                 currentDoF_kneeJoint_gt_part1; ...
                                 currentDoF_kneeJoint_gt_part2])]);


        % 4) Extra --------------------------------------------------------

        [h, p] = ttest2(currentDoF_kneeJoint_estgt_diff1, currentDoF_kneeJoint_estgt_diff2, 'Vartype','unequal');
        pvalues_rmse(idx_dof, :) = [h, p];

        [p, h] = ranksum(abs(currentDoF_kneeJoint_estgt_diff1), abs(currentDoF_kneeJoint_estgt_diff2));
        pvalues_mad(idx_dof, :) = [h, p];

end


% 4) Swap row-4 and row-6, from r-zyx to r-xyz
tmp = metric_mat_part1(4, :);
metric_mat_part1(4, :) = metric_mat_part1(6, :);
metric_mat_part1(6, :) = tmp;

% 5) Create a table. You can just open the table and copy paste to excel 
% or latex table generator.
metric_table_part1 = array2table( metric_mat_part1,  'VariableNames', ...
                                  {'rmse_mean', 'rmse_std', 'mad_q2', 'mad_uw', 'corr'});

% Repeat 4, for part 2
tmp = metric_mat_part2(4, :);
metric_mat_part2(4, :) = metric_mat_part2(6, :);
metric_mat_part2(6, :) = tmp;

% Repeat 5, for part 2
metric_table_part2 = array2table( metric_mat_part2,  'VariableNames', ...
                                  {'rmse_mean', 'rmse_std', 'mad_q2', 'mad_uw', 'corr'});

% Repeat 4 for pvalues_rmse
tmp = pvalues_rmse(4, :);
pvalues_rmse(4, :) = pvalues_rmse(6, :);
pvalues_rmse(6, :) = tmp;

% Repeat 4 for pvalues_mad
tmp = pvalues_mad(4, :);
pvalues_mad(4, :) = pvalues_mad(6, :);
pvalues_mad(6, :) = tmp;