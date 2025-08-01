%% HELOOOOOO
% - Before running this script, you must already generated .mat files:
%   1) all_TsReg_sX_mXX.mat
%      - Contains transformations of the bone (tibia) relative to the ref 
%        (global) both GT and est.
%      - Generated from main3_registrationWithTime_4clean.m
%   2) all_kneeJoint6DOFs_sX_mXX.mat
%      - Contains transformations of the bones (femur and tibia) relative to 
%        the ref (global), and the 6Dof of knee joint kinematic. 
%      - Generated from main4_kinematicEstimation3_TFJointRUMC.m
%   3) cycle_timestamp_sX_mXX.csv
%      - Contains the timestamps of the cycle. 
%       - Generated from extra_detectFErotCycle.m
% - The purpose of this script is to quantify and evaluate the 3d pose 
%   estimation against the ground truth.
% - main5b_registrationEval2 will show multiple registration errors from
%   different method at once

clc; clear; close all;

% [EDIT] directory to the project
path_root    = 'D:\Documents\BELANDA\PhD Thesis\Code\MATLAB\amode_navigation_experiment\experiment_b';

% [EDIT] Change the data you are using accordingly. 
% -----> dir_depthdata is created by main1_processDepthData.m
% -----> dir_Tdata is created by main3_registrationWithTime.m
dirs_Tdata = { fullfile('depthdata_s4_m04_20250708-172830', 'Tdata_s4_m04_20250724-100122'), ...   % with-nav
               fullfile('depthdata_s3_m02_20250722-114731', 'Tdata_s3_m02_20250724-030951'), ...   % no-nav, manual
               fullfile('depthdata_s3_m02_20250722-174503', 'Tdata_s3_m02_20250724-032542')};      % no-nav, auto, 2x noise

% [EDIT] select which cycle you want to show
cycle_select  = [2,5];

% [EDIT]
% do you want to split the data into cycles and average them?
is_split  = false;
% if true, do you want to show the std as a shade or washout-lines?
is_shaded = false;

% [EDIT] select the color scheme
color_scheme = 2;

% [EDIT] for saving the resulting mat file
is_saveMat = false;



%% INITIALIZE FIGURE OBJECTS AND EVERYTHING RELATED TO PLOTS

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

% bar chart
fig1 = figure('Name', 'Registration Error', 'Position', [50 50 1100 500]);
% create a tab group
tg1 = uitabgroup(fig1);
% create tab and axis
tabtitle = {'us2regbone', 'regbone2gtbone'};
ylims    = [3, 5];
ax1 = {};
for idx_metric=1:2
    % make a new tab
    tb1 = uitab(tg1, 'Title', tabtitle{idx_metric});

    % set the properties of the axes
    ax_tmp = axes(tb1);
    hold(ax_tmp, 'on');
    grid(ax_tmp, 'on');
    axis(ax_tmp, 'tight');
    title(ax_tmp, 'Mean Registration Error Through Time');
    xlabel(ax_tmp, 'Time (s)');
    ylabel(ax_tmp, 'Error (mm)');
    ylim(ax_tmp, [0 ylims(idx_metric)]);
    ax_tmp.FontSize = 12;

    % store the axes
    ax1{idx_metric} = ax_tmp;
end

% bar chart
fig2 = figure('Name', 'Registration Error');
ax2 = axes(fig2);
hold(ax2, 'on');
grid(ax2, 'on');
axis(ax2, 'tight');
title(ax2, 'Registration error');
xlabel(ax2, 'Distance Metric');
ylabel(ax2, 'RMSE (mm)');
ax2.FontSize = 12;



%% INITIALIZE PATHS AND LOADING SOME CONFIGURATION

% declare some of the important paths
path_function     = fullfile(path_root, 'functions');
path_outputs      = fullfile(path_root, 'outputs');
path_data         = fullfile(path_root, 'data');

% Generate path to function directory
addpath(genpath(path_function));


%% INITIALIZE SOME DATA
% - Why do we need to do this in advance before we plot our data? Because we
%   are loading different methods (from different file) that have different
%   cycle duration pattern. 
% - We can't do something like load data -> plot -> load another -> plot,
%   because the data is not matching.
% - To match all of the data, we must load all the data first and compute the
%   length and make a strategy to cut the data.

n_dirs_Tdata = length(dirs_Tdata);
all_kneeJoint6DOFs_tables = cell(1, n_dirs_Tdata);
all_TsReg_tables          = cell(1, n_dirs_Tdata);
all_cycle_timestamps      = cell(1, n_dirs_Tdata);

% loop trough all method we have (different .mat file)
for idx_dir=1:n_dirs_Tdata

    tmp_str  = strsplit(dirs_Tdata{idx_dir}, {'_', '\'});
    sess_str = tmp_str{2};
    meas_str = tmp_str{3};

    % 1.a) Load all the data related to knee joint 6 dof estimation.
    % -> This mat file is generated from main4_kinematicEstimation.m
    mat_filename = sprintf('all_kneeJoint6DOFs_%s_%s.mat', sess_str, meas_str);
    mat_fullpath = fullfile(path_outputs, 'output_allest', dirs_Tdata{idx_dir}, mat_filename);
    load(mat_fullpath);
    % 1.b) Store the table to our cells
    all_kneeJoint6DOFs_tables{idx_dir} = all_kneeJoint6DOFs_table;
    
    
    % 2.a) Load all the data related to bone pose estimation performance
    %   -> This mat file is generated from main3_registrationWithTime.m
    mat_filename = sprintf('all_TsReg_%s_%s.mat', sess_str, meas_str);
    mat_fullpath = fullfile(path_outputs, 'output_allest', dirs_Tdata{idx_dir}, mat_filename);
    load(mat_fullpath);
    % 2.b) Store the table to our cells
    all_TsReg_tables{idx_dir} = all_TsReg_table;
    
    % 3) Load the cycle timestamp data. 
    % -> This data indicates where the cycle motion starts, what was been 
    %    done in the experiment.
    % -> this is obtained from extra_detectFErotCycle.m
    csv_filename = sprintf('cycle_timestamp_%s_%s.csv', sess_str, meas_str);
    csv_fullpath = fullfile(path_outputs, csv_filename);
    cycle_timestamp = readmatrix(csv_fullpath);
    % 3.b) Store the cycle to our cells
    all_cycle_timestamps{idx_dir} = cycle_timestamp;

end


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


%% Figure 1

style = [];
style(1).alphaVal  = 1.0;
style(1).lineWidth = 2.0;
style(1).lineStyle = '-';
style(2).alphaVal  = 0.6;
style(2).lineWidth = 1.0;
style(2).lineStyle = '--';
style(3).alphaVal  = 0.8;
style(3).lineWidth = 1.0;
style(3).lineStyle = ':';

% There are two different registration metric, us2regbone and regbone2gtbone
% so i will loop this process two times. The only difference is just:
% 1) the time i retrieve the registration error data
% 2) plotting, i need to put it to different axes object
for idx_metric=1:2

    % Loop trough all method we have (different .mat file)
    for idx_dir=1:n_dirs_Tdata

        % 1) Get the current data for cycle_timestamps
        cycle_timestamps         = all_cycle_timestamps{idx_dir};
        cycle_idxoriginal_select = [cycle_timestamps(cycle_select(1)), cycle_timestamps(cycle_select(2))];
    
        % 2) Get the current data for knee joint data (and registration data)
        all_kneeJoint6DOFs_table = all_kneeJoint6DOFs_tables{idx_dir};
        all_TsReg_table          = all_TsReg_tables{idx_dir};
    
        % -> Delete problematic rows . Here, I am only using 
        %    all_kneeJoint6DOFs_table to get the is_invalid column. 
        % -> I am assuming that all_kneeJoint6DOFs_table and all_TsReg_table 
        %    are synchronized (as it should be). 
        % -> I tried to be very careful when when i generated those data.
        idcs_problematicRows = find(all_kneeJoint6DOFs_table.is_invalid);
        all_kneeJoint6DOFs_table(idcs_problematicRows, :) = [];
        all_TsReg_table(idcs_problematicRows, :) = [];
    
        % 3) Grab the valid index (same assumption, i assume the syncs) 
        timestamp_idcs_valid = all_kneeJoint6DOFs_table.Timestamp_idx;
        timestamp_ms_valid = all_kneeJoint6DOFs_table.Timestamp_ms;  
    
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
    
        % 5) Grab the knee joint 6dof and convert it into matrix. 
        % -> Here i only take the ground truth. I just want to use this data 
        %    for reference, to find a pattern within registration performance
        kneeJoint6DOFs_gt       = cell2mat(all_kneeJoint6DOFs_table.kneeJoint6DOFs_gt);
        % -> Get the registration errors (it is already a vector, no need convert)
        if(idx_metric==1)
            errors_rmse_us2regbone  = all_TsReg_table.errors_rmse_us2regbone;
        else
            errors_rmse_us2regbone  = all_TsReg_table.errors_rmse_regbone2gtbone;
        end
    
        % Note ------------------------------------------------------------
        % -> Here i don't need to get the kneeJoint est data)
        % -> And, i also don't need to load all the dof, i only focus on
        %    flexion-extension)

        % 2) Get the current gt dof
        selDoF_kneeJoint_gt    = kneeJoint6DOFs_gt(cycle_idxvalid_select(1):cycle_idxvalid_select(2), 4);

        % 3) Get the regisration data
        selErr_rmse_us2regbone = errors_rmse_us2regbone(cycle_idxvalid_select(1):cycle_idxvalid_select(2));
        % -> and the selected timestamp
        selTime_ms_valid       = timestamp_ms_valid(cycle_idxvalid_select(1):cycle_idxvalid_select(2));
    
        % 4) Stretch the value, here i am using imresize from image processing
        %    tool box, because why not? they already have built-in
        selDoF_kneeJoint_gt_valuestretched    = imresize(selDoF_kneeJoint_gt, [max_timestamps 1], 'bicubic');
        selErr_rmse_us2regbone_valuestretched = imresize(selErr_rmse_us2regbone, [max_timestamps 1], 'bicubic');
        selTime_ms_valid_valuestretched       = imresize(selTime_ms_valid, [max_timestamps 1], 'bicubic');
    
        % 5) Plot the data
        % -> Okay here is my mistake. I am building this particular
        %    function that plot the registration error trough time with
        %    color coded that corresponds to the angle of flex-ext.
        % -> In order to create colorful line, I created tons of lines to
        %    connect each two consecutive data points.
        % -> By doing this, we can't easily set the x-axis to the timestamp
        %    that we select and process before
        % -> Because of this, there is an unnecessary process after the loop
        plotColoredTrace( selDoF_kneeJoint_gt_valuestretched', selErr_rmse_us2regbone_valuestretched', ax1{idx_metric}, ...
                          style(idx_dir).alphaVal, style(idx_dir).lineWidth, style(idx_dir).lineStyle);
    
    end

    % (Here, I should add legend. But I can't. The way i form the colored
    % line is by creating hundreds of line with different color. Legend 
    % will only show some of them. My solution is, leave it like this, you
    % can later edit the figure yourself with Adobe Illustrator)

    % Adjust the x-axis of the plot. 
    % -> I should have done this when i plot the lines in the plotColoredTrace 
    %    function. But it was quite complicated because i plot it 2-lines at a 
    %    time (to support colour). 
    % -> So better i arrange the x-axis here instead of there. 
    % -> Here, I am create my own x-tick according the timestamp_ms values
    tmp_end    = round(selTime_ms_valid_valuestretched(end) - selTime_ms_valid_valuestretched(1));
    tmp_tick   = round(tmp_end / (length(ax1{idx_metric}.XTick)+1));
    
    % Create the xtick
    tmp_xticks = 0:tmp_tick:tmp_end;       % [ms]
    tmp_xticks = round(tmp_xticks/ 1000);  % [s]
    
    % Replace the xtick string/label
    ax1{idx_metric}.XTickLabel = tmp_xticks(2:end);

% bye-bye, it works, i don't care.
end

%% FIGURE 2

% prepare the variable for the dabarplot
all_RMSE_value  = [];
all_RMSE_grpidx = [];
all_RMSE_grpstr = {};
xtlabels_str = {'us2estbone', 'estbone2gtbone'};

% Prepare the variables for the statistic table.
% rows will be organized for the metrics, see xtlabels_str.
% cols will be organized as: (rmse, std) times the number of n_dirs_Tdata.
%   -> Be aware the orders. Please follow all_RMSE_grpstr.
metric_mat = zeros(length(xtlabels_str), 2*n_dirs_Tdata);

% [EDIT] The names of the methods.
%     -> this value should be adjusted to how many method you use
all_RMSE_grpstr = {'No Navigation, manual', 'No Navigation, auto', 'With Navigation'};

% [EDIT] To arrange the order, it is better from worst to best
%     -> same as above, need to be the same number of methods you use
bargroup_showorder = [3, 2, 1];

% Loop trough all method we have (different .mat file)
for idx_dir=1:n_dirs_Tdata

    % 1) Get the current data for cycle_timestamps
    cycle_timestamps         = all_cycle_timestamps{idx_dir};
    cycle_idxoriginal_select = [cycle_timestamps(cycle_select(1)), cycle_timestamps(cycle_select(2))];

    % 2) Get the current data for registration data (similar step to all_kneeJoint6DOFs_table)
    all_TsReg_table          = all_TsReg_tables{idx_dir};
    % -> delete problematic rows
    idcs_problematicRows = find(all_TsReg_table.is_invalid);
    all_TsReg_table(idcs_problematicRows, :) = [];

    % 3) Get the valid timestamp.
    % -> As i mentioned in the first part of the script (figure 1), here 
    %    I am assuming that all_kneeJoint6DOFs_table and all_TsReg_table 
    %    are synchronized. 
    timestamp_idcs_valid = all_TsReg_table.Timestamp_idx;
    timestamp_ms_valid = all_TsReg_table.Timestamp_ms;

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

    % ---------------------------------------------------------------------
    % (Here is a bit different than the usual steps. Because here we want
    % to plot a barplot with dabarplot function) 

    % 1) Arrange the values
    RMSE_value  = [ all_TsReg_table.errors_rmse_us2regbone(cycle_idxvalid_select(1):cycle_idxvalid_select(2)), ...
                    all_TsReg_table.errors_rmse_regbone2gtbone(cycle_idxvalid_select(1):cycle_idxvalid_select(2)) ];
    RMSE_grpidx = bargroup_showorder(idx_dir) .* ones(1, size(RMSE_value, 1));

    % 2) Store the values
    all_RMSE_value = [all_RMSE_value; RMSE_value];
    all_RMSE_grpidx = [all_RMSE_grpidx, RMSE_grpidx];


    % Compute the mean and std for metric mat
    tmp_rmse = mean(RMSE_value, 1);
    tmp_std  = std(RMSE_value, 1);
    metric_mat( :, bargroup_showorder(idx_dir)*2 -1) = tmp_rmse';
    metric_mat( :, bargroup_showorder(idx_dir)*2   ) = tmp_std';

end

% Activate the corresponding axis (daboxplot, for some reason, does not 
% have argument to select particular axes)
axes(ax2);

% Draw the bar plot
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

% Create a statistic/metric table. You can just open the table and copy 
% paste to excel or latex table generator.
metric_table = array2table( metric_mat,  'VariableNames', ...
                            { 'method1_mean', 'method1_std', ...
                              'method2_mean', 'method2_std', ...
                              'method3_mean', 'method3_std'});
% I add this variable, just in case if dennis in the future don't know what
% the order of method1, method2, method3
metric_table_methodstring = all_RMSE_grpstr;