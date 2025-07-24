clc; clear; close all;

% [EDIT] directory to the project
path_root    = 'D:\DennisChristie\NavGuidedAmodeUS-BoneKinematicProcessing';

% [EDIT] Change the data you are using accordingly. 
% -----> dir_depthdata is created by main1_processDepthData.m
% -----> dir_Tdata is created by main3_registrationWithTime.m
dirs_Tdata = { fullfile('depthdata_s4_m04_20250708-172830', 'Tdata_s4_m04_20250724-100122'), ...   % with-nav
               fullfile('depthdata_s3_m02_20250722-114731', 'Tdata_s3_m02_20250724-030951'), ...   % no-nav, manual
               fullfile('depthdata_s3_m02_20250722-174503', 'Tdata_s3_m02_20250724-032542')};      % no-nav, auto, 2x noise


% [EDIT] Color scheme.
color_scheme = 2;


%% INITIALIZE FIGURE OBJECTS AND EVERYTHING RELATED TO PLOTS

% set the ylim so that the axes become consistent
ax_ylim = repmat([-20, 20], 6, 1);

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


% First figure is to show the joint kinematic with all of the cycle parts
fig1 = figure('Name', 'Joint Kinematic: Relative to Ground Truth', 'Position', [50 50 700 500]);
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
    xlabel(ax1(idx_dof), 'Timestamp');
    if(idx_dof<=3)
        ylabel(ax1(idx_dof), 'mm');
    else
        ylabel(ax1(idx_dof), 'deg');
    end
    ax1(idx_dof).FontSize = 12;
    ylim( ax1(idx_dof), ax_ylim(idx_dof,:) );
end


%% INITIALIZE PATHS AND LOADING SOME CONFIGURATION

% declare some of the important paths
path_function     = fullfile(path_root, 'functions');
path_outputs      = fullfile(path_root, 'outputs');
path_data         = fullfile(path_root, 'data');

% Generate path to function directory
addpath(genpath(path_function));


%% INITIALIZE SOME DATA

n_dirs_Tdata = length(dirs_Tdata);
all_kneeJoint6DOFs_tables = cell(1, n_dirs_Tdata);
all_cycle_timestamps       =  cell(1, n_dirs_Tdata);

for idx_dir=1:n_dirs_Tdata
    % 1) Load all the data related to knee joint 6 dof estimation.
    % -> This mat file is generated from main4_kinematicEstimation.m
    tmp_str  = strsplit(dirs_Tdata{idx_dir}, {'_', '\'});
    sess_str = tmp_str{2};
    meas_str = tmp_str{3};
    mat_filename = sprintf('all_kneeJoint6DOFs_%s_%s.mat', sess_str, meas_str);
    mat_fullpath = fullfile(path_outputs, 'output_allest', dirs_Tdata{idx_dir}, mat_filename);
    load(mat_fullpath);
    % 1.b) Store the table to our cells
    all_kneeJoint6DOFs_tables{idx_dir} = all_kneeJoint6DOFs_table;
    
    % 2) Load the cycle timestamp data. 
    % -> This data indicates where the cycle motion starts, what was been 
    %    done in the experiment.
    % -> this is obtained from extra_detectFErotCycle.m
    csv_filename = sprintf('cycle_timestamp_%s_%s.csv', sess_str, meas_str);
    csv_fullpath = fullfile(path_outputs, csv_filename);
    cycle_timestamps = readmatrix(csv_fullpath);
    % 2.b) Store the cycle to our cells
    all_cycle_timestamps{idx_dir} = cycle_timestamps;
end



%% EXTRA

% select cycle
cycle_select             = [2,5];
% to store the max timestamps
max_timestamps = 0;

% check which data has the longest timestamps for all cycle
for idx_dir=1:n_dirs_Tdata

    % get the current data (timestamps)
    cycle_timestamps = all_cycle_timestamps{idx_dir};

    % get the length of the selected cycles
    cycle_idxoriginal_select   = ( cycle_timestamps(cycle_select(1)):cycle_timestamps(cycle_select(2)) );
    n_timestamps_allcycle      = length(cycle_idxoriginal_select);

    % store the number if the current data  is longer than before
    if(n_timestamps_allcycle>max_timestamps)
        max_timestamps = n_timestamps_allcycle;
    end
end




%%

for idx_dir=1:n_dirs_Tdata

    % get the current data
    all_kneeJoint6DOFs_table = all_kneeJoint6DOFs_tables{idx_dir};
    cycle_timestamps = all_cycle_timestamps{idx_dir};

    % delete problematic rows
    idcs_problematicRows = find(all_kneeJoint6DOFs_table.is_invalid);
    all_kneeJoint6DOFs_table(idcs_problematicRows, :) = [];

    % Grab the valid index
    timestamp_idcs_valid = all_kneeJoint6DOFs_table.Timestamp_idx;
    timestamp_ms_valid = all_kneeJoint6DOFs_table.Timestamp_ms;
    
    % Grab the knee joint 6dof and convert it into matrix
    kneeJoint6DOFs_est = cell2mat(all_kneeJoint6DOFs_table.kneeJoint6DOFs_est);
    kneeJoint6DOFs_gt  = cell2mat(all_kneeJoint6DOFs_table.kneeJoint6DOFs_gt);

    % get the cycle
    cycle_idxoriginal_select = [cycle_timestamps(cycle_select(1)), cycle_timestamps(cycle_select(2))];

    % The indices that are shown in cycle_timestamp is the indices of the table
    % original measurement. We make a lot of filtering process already
    % (removing some invalid data). So we need to get the actual timestamp.
    tmp1 = find(timestamp_idcs_valid==cycle_idxoriginal_select(1));
    tmp2 = find(timestamp_idcs_valid==cycle_idxoriginal_select(2));
    cycle_idxvalid_select = [tmp1, tmp2];

    % loop for all dofs
    for idx_dof=1:6
    
        % get the current est dof
        currentDoF_kneeJoint_est = kneeJoint6DOFs_est(cycle_idxvalid_select(1):cycle_idxvalid_select(2), idx_dof);
        currentDoF_kneeJoint_est = smoothdata(currentDoF_kneeJoint_est, 'rlowess', 20);

        % get the current gt dof
        currentDoF_kneeJoint_gt = kneeJoint6DOFs_gt(cycle_idxvalid_select(1):cycle_idxvalid_select(2), idx_dof);
    
        % calculate the difference
        currentDoF_kneeJoint_estgtdiff   = currentDoF_kneeJoint_est - currentDoF_kneeJoint_gt;

        % stretch the value, here i am using imresize from image processing
        % tool box, because why not? they already have built-in
        currentDoF_kneeJoint_est_valuestretched = imresize(currentDoF_kneeJoint_est, [max_timestamps 1], 'bicubic');
        currentDoF_kneeJoint_estgtdiff_valuestretched = imresize(currentDoF_kneeJoint_estgtdiff, [max_timestamps 1], 'bicubic');
    
        % plot the knee joint for the current dof
        plot(ax1(idx_dof), currentDoF_kneeJoint_estgtdiff_valuestretched, '-', 'Color', c(idx_dir+1, :), 'LineWidth', 3);
    
        % compute the quantitative result
        result_corr  = corr(currentDoF_kneeJoint_est, currentDoF_kneeJoint_gt);
        result_rmse  = sqrt( mean( (currentDoF_kneeJoint_estgtdiff).^2 ) );
        result_std   = std(currentDoF_kneeJoint_estgtdiff);
        [~, ~, ~, result_q2, ~, result_uw, ~] = computeBoxplotStats(abs(currentDoF_kneeJoint_estgtdiff));
    
    end

end

%% COSMETIC

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
                'Distraction(+) - Compression(-)', ...
                'Medial(+) - Lateral(-)', ...
                'Flexion(+) - Extension(-)', ...
                'External(+) - Internal(-)', ...
                'Abduction(+) - Adduction(-)', };
end

for idx_dof = 1:6
    title(ax1(idx_dof), ax_title{idx_dof});
    plot(ax1(idx_dof), zeros(length(currentDoF_kneeJoint_estgtdiff_valuestretched), 1), '-', 'Color', 'b', 'LineWidth', 3);
    legend(ax1(idx_dof), {'with-nav', 'no-nav, manual', 'no-nav, 2x-noise', 'Ground truth'});
end



















