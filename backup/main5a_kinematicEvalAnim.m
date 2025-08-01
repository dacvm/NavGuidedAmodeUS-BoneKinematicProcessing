%% HELOOOOOO
% - Before running this script, you must already generated .mat file called
%   "all_kneeJoint6DOFs_sX_mXX.mat". This .mat file contains transformations 
%   of the bones (femur and tibia) relative to the ref (global) both GT 
%   and est.
% - This script will produce 6 plots, each plot for each individual degree
%   of freedom of knee joint kinematic estimation against its ground truth
% - Note, that the calculation of joint knee kinematic is defined 
%   in main4_kinematicEstimation.m
% - There are two options to show the plot: all continuous cycle shown at
%   once, or cycles are chopped and each chop shown overalayed to each
%   other together with their mean.
% - If you want to generate the error relative to ground truth, cehck
%   main5a_kinematicEvalAll.m

clc; clear; close all;

% [EDIT] directory to the project
path_root    = 'D:\Documents\BELANDA\PhD Thesis\Code\MATLAB\amode_navigation_experiment\experiment_b';

% [EDIT] Change the data you are using accordingly. 
% -----> dir_depthdata is created by main1_processDepthData.m
% -----> dir_Tdata is created by main3_registrationWithTime.m
dirs_Tdata = { fullfile('depthdata_s4_m04_20250708-172830', 'Tdata_s4_m04_20250724-100122'), ...   % with-nav
               fullfile('depthdata_s3_m02_20250722-114731', 'Tdata_s3_m02_20250724-030951'), ...   % no-nav, manual
               fullfile('depthdata_s3_m02_20250722-174503', 'Tdata_s3_m02_20250724-032542')};      % no-nav, auto, 2x noise

% [EDIT] select which data you want to show
idx_dir_Tdata = 1;

% [EDIT] select which cycle you want to show
cycle_select  = [2,5];

% [EDIT] select the color scheme
color_scheme = 2;

% [EDIT] for saving the resulting mat file
is_saveMat = false;


%% INITIALIZE PATHS AND LOADING SOME CONFIGURATION

% declare some of the important paths
path_function     = fullfile(path_root, 'functions');
path_outputs      = fullfile(path_root, 'outputs');
path_data         = fullfile(path_root, 'data');

% declare some of important path of the data
path_bonestl  = fullfile(path_root, "data", "ct", "bone");

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

% 3) Load all the data related to CT
run('extra_structCTdata.m');

% 4) Create base rotation to trasnform qualisys base vector to matlab
% -- Qualisys has y direction as up, MATLAB has z direction as up
R_tmp = eul2rotm([0 0 pi/2], "ZYX");
t_tmp = [0 0 0]';
baseRotation_Qualisys2Matlab = [R_tmp, t_tmp; 0 0 0 1];

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
                'Distraction(+) - Compression(-)', ...
                'Medial(+) - Lateral(-)', ...
                'Flexion(+) - Extension(-)', ...
                'External(+) - Internal(-)', ...
                'Abduction(+) - Adduction(-)', };
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
fig1 = figure('Name', 'Joint Kinematic: All Cycle Parts', 'WindowState', 'maximized');
t1 = tiledlayout(fig1, 4, 3, ...
     'TileSpacing', 'compact', ...   % tighten spacing if you like
     'Padding',     'compact');      % remove outer margins

% Pre-allocate an array of axes handles
ax1 = gobjects(3,1);

% Populate each tile and store its axes handle
for idx_view = 1:3
    ax1(idx_view) = nexttile(t1, idx_view, [2 1]);
    hold(ax1(idx_view), 'on');
    grid(ax1(idx_view), 'on');
    grid(ax1(idx_view), 'minor');
    axis(ax1(idx_view), 'equal');
    axis(ax1(idx_view), 'tight');
    xlabel(ax1(idx_view), '← Anterior | X (mm) | Posterior →');
    ylabel(ax1(idx_view), '← Lateral | Y (mm) | Medial →');
    zlabel(ax1(idx_view), '← Distal | Z (mm) | Proximal →');
    xlim(ax1(idx_view), [-60 60]);
    ylim(ax1(idx_view), [-60 60]);
    zlim(ax1(idx_view), [-60 60]);
    ax1(idx_view).FontSize = 12;
end

title(ax1(1), 'Frontal Plane (from front)');
view(ax1(1), 90, 0);
title(ax1(2), 'Sagital Plane (from lateral)');
view(ax1(2), 0, 0);
title(ax1(3), 'Transverse Plane (from top)');
view(ax1(3), 0, 90);


% Pre-allocate an array of axes handles
ax2 = gobjects(6,1);

% Populate each tile and store its axes handle
for idx_dof = 1:6
    ax2(idx_dof) = nexttile(t1, idx_dof + (length(ax1)*2) ) ;
    hold(ax2(idx_dof), 'on');
    grid(ax2(idx_dof), 'on');
    axis(ax2(idx_dof), 'tight');
    title(ax2(idx_dof), ax_title{idx_dof});
    xlabel(ax2(idx_dof), 'Timestamp');
    if(idx_dof<=3)
        ylabel(ax2(idx_dof), 'mm');
    else
        ylabel(ax2(idx_dof), 'deg');
    end
    ax2(idx_dof).FontSize = 12;
    % ylim( ax1(idx_dof), ax_ylim(idx_dof,:) );
end



%%

% delete problematic rows
idcs_problematicRows = find(all_kneeJoint6DOFs_table.is_invalid);
all_kneeJoint6DOFs_table(idcs_problematicRows, :) = [];

% Grab the valid index
timestamp_idcs_valid = all_kneeJoint6DOFs_table.Timestamp_idx;
timestamp_ms_valid = all_kneeJoint6DOFs_table.Timestamp_ms;

% Grab the knee joint 6dof and convert it into matrix
kneeJoint6DOFs_est = cell2mat(all_kneeJoint6DOFs_table.kneeJoint6DOFs_est);
kneeJoint6DOFs_gt  = cell2mat(all_kneeJoint6DOFs_table.kneeJoint6DOFs_gt);

% select cycle
cycle_idxoriginal_select = [cycle_timestamp(cycle_select(1)), cycle_timestamp(cycle_select(2))];

% The indices that are shown in cycle_timestamp is the indices of the table
% original measurement. We make a lot of filtering process already
% (removing some invalid data). So we need to get the actual timestamp.
tmp1 = find(timestamp_idcs_valid==cycle_idxoriginal_select(1));
tmp2 = find(timestamp_idcs_valid==cycle_idxoriginal_select(2));
cycle_idxvalid_select = [tmp1, tmp2];

% 1) First we need to show he whole plot first ============================

% variable to store the smoothed data
kneeJoint6DOFs_estSmooth = zeros(size(kneeJoint6DOFs_est));

% loop for all dofs, let's display them one by one
for idx_dof=1:6

    % get the current dof
    currentDoF_kneeJoint_est = kneeJoint6DOFs_est(cycle_idxvalid_select(1):cycle_idxvalid_select(2), idx_dof);
    currentDoF_kneeJoint_est = smoothdata(currentDoF_kneeJoint_est, 'rlowess', 20);
    currentDoF_kneeJoint_gt = kneeJoint6DOFs_gt(cycle_idxvalid_select(1):cycle_idxvalid_select(2), idx_dof);

    % calculate the difference
    currentDoF_kneeJoint_estgt_diff   = currentDoF_kneeJoint_est - currentDoF_kneeJoint_gt;

    % plot the knee joint for the current dof
    plot(ax2(idx_dof), cycle_idxvalid_select(1):cycle_idxvalid_select(2), currentDoF_kneeJoint_est, '-', 'Color', 'r', 'LineWidth', 2);
    plot(ax2(idx_dof), cycle_idxvalid_select(1):cycle_idxvalid_select(2), currentDoF_kneeJoint_gt, '-', 'Color', 'b', 'LineWidth', 2);
    % legend(ax2(idx_dof), {'Estimation', 'Ground truth'});

    % compute the quantitative result
    result_corr  = corr(currentDoF_kneeJoint_est, currentDoF_kneeJoint_gt);
    result_rmse  = sqrt( mean( (currentDoF_kneeJoint_estgt_diff).^2 ) );
    result_std   = std(currentDoF_kneeJoint_estgt_diff);
    [~, ~, ~, result_q2, ~, result_uw, ~] = computeBoxplotStats(abs(currentDoF_kneeJoint_estgt_diff));

    str = sprintf('corr   = %.2f\nrmse = %.2f, %.2f\nmad  = %.2f%s%.2f', result_corr, result_rmse, result_std, result_q2, char(9652), result_uw);
    text( ax2(idx_dof), ...
          0.05, 0.95, ... 
          str, ...
          'Units',            'normalized', ...
          'VerticalAlignment','top', ...
          'FontSize',         10, ...
          'BackgroundColor',  'w', ...
          'EdgeColor',        'k', ...
          'Margin',           5 );

    % store the value
    kneeJoint6DOFs_estSmooth(cycle_idxvalid_select(1):cycle_idxvalid_select(2), idx_dof) = currentDoF_kneeJoint_est;

end


% 2) Second, Animation ====================================================
% loop for all the time stamps
for idx_t = cycle_idxvalid_select(1):cycle_idxvalid_select(2)

    % 2.0) Delete previous objects in the plot in the provided axes
    for i=1:length(ax1)
        delete(findobj(ax1(i), 'Tag', '3d_bone_gt'));
        delete(findobj(ax1(i), 'Tag', 'cs_bone_gt'));
        delete(findobj(ax1(i), 'Tag', 'text_t_gt'));
        delete(findobj(ax1(i), 'Tag', 'text_r_gt'));
        delete(findobj(ax1(i), 'Tag', 'line_t_gt'));
    end

    for i=1:length(ax2)
        delete(findobj(ax2(i), 'Tag', 'current_time'));
        delete(findobj(ax2(i), 'Tag', 'current_gt'));
        delete(findobj(ax2(i), 'Tag', 'current_est'));
    end


    % 2.1) Femur ----------------------------------------------------------

    % Get the bone coordinate frame relative to CT scan
    T_bone_ct    = allBone_CT(1).T;

    % Get the femur coordinate frame relative to itself (well, identity)
    T_bone_femur = eye(4);   

    % Transform the point cloud
    current_boneCTpoints_GT = baseRotation_Qualisys2Matlab * T_bone_femur * inv(T_bone_ct) * ...
                              [allBone_CT(1).stl.Points'; ones(1, size(allBone_CT(1).stl.Points, 1))];
    % Construct the mesh
    current_boneCTtri_GT    = triangulation( allBone_CT(1).stl.ConnectivityList, ...
                                             current_boneCTpoints_GT(1:3,:)');

   % Display the bone and the CS on all of the 3-view axes
    for i=1:3
        % display the surface
        trisurf(current_boneCTtri_GT, 'FaceColor', '#bdc3c7', 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'EdgeAlpha', 0.02, 'Parent', ax1(i), 'Tag', '3d_bone_gt');
        % display the coordinate frame
        tmp = baseRotation_Qualisys2Matlab * T_bone_femur;
        display_axis_v2(ax1(i), tmp(1:3, 4), tmp(1:3, 1:3), 15, '', 'Tag', 'cs_bone_gt');
    end


    % 2.2) Tibia GT -------------------------------------------------------

    % Get the bone coordinate frame relative to CT scan
    T_bone_ct    = allBone_CT(2).T;

    % Get the tibia est coordinate frame relative to the Femur
    T_femurGT_ref = all_kneeJoint6DOFs_table.Ts_femurGT_ref{idx_t};
    T_tibiaGT_ref = all_kneeJoint6DOFs_table.Ts_tibiaGT_ref{idx_t};
    T_bone_femur  = inv(T_femurGT_ref) * T_tibiaGT_ref;

    % Transform the point cloud
    current_boneCTpoints_GT = baseRotation_Qualisys2Matlab * T_bone_femur * inv(T_bone_ct) * ...
                              [allBone_CT(2).stl.Points'; ones(1, size(allBone_CT(2).stl.Points, 1))];
    % Construct the mesh
    current_boneCTtri_GT    = triangulation( allBone_CT(2).stl.ConnectivityList, ...
                                             current_boneCTpoints_GT(1:3,:)');

   % Display the bone and the CS on all of the 3-view axes
    for i=1:3
        % display the surface
        trisurf(current_boneCTtri_GT, 'FaceColor', '#bdc3c7', 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'EdgeAlpha', 0.02, 'Parent', ax1(i), 'Tag', '3d_bone_gt');
        % display the coordinate frame
        tmp = baseRotation_Qualisys2Matlab * T_bone_femur;
        display_axis_v2(ax1(i), tmp(1:3, 4), tmp(1:3, 1:3), 15, '', 'Mode', 'thin', 'Tag', 'cs_bone_gt');
    end


    % 2.3) Tibia Est ------------------------------------------------------

    % Get the bone coordinate frame relative to CT scan
    T_bone_ct    = allBone_CT(2).T;

    % Get the tibia est coordinate frame relative to the Femur
    T_femurGT_ref  = all_kneeJoint6DOFs_table.Ts_femurGT_ref{idx_t};
    T_tibiaEst_ref = all_kneeJoint6DOFs_table.Ts_tibiaEst_ref{idx_t};
    T_bone_femur  = inv(T_femurGT_ref) * T_tibiaEst_ref;

    % Transform the point cloud
    current_boneCTpoints_GT = baseRotation_Qualisys2Matlab * T_bone_femur * inv(T_bone_ct) * ...
                              [allBone_CT(2).stl.Points'; ones(1, size(allBone_CT(2).stl.Points, 1))];
    % Construct the mesh
    current_boneCTtri_GT    = triangulation( allBone_CT(2).stl.ConnectivityList, ...
                                             current_boneCTpoints_GT(1:3,:)');

   % Display the bone and the CS on all of the 3-view axes
    for i=1:3
        % display the surface
        trisurf(current_boneCTtri_GT, 'FaceColor', '#2ecc71', 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'EdgeAlpha', 0.08, 'Parent', ax1(i), 'Tag', '3d_bone_gt');
        % display the coordinate frame
        tmp = baseRotation_Qualisys2Matlab * T_bone_femur;
        display_axis_v2(ax1(i), tmp(1:3, 4), tmp(1:3, 1:3), 15, '', 'Tag', 'cs_bone_gt');
    end

    % Test ----------------------------------------------------------------

    tmp = baseRotation_Qualisys2Matlab * T_bone_femur;
    t_tibia = tmp(1:3, 4);
    R_tibia = tmp(1:3, 1:3);
    r_tibia = rad2deg(rotm2eul(R_tibia));
    r_tibia = [r_tibia(1), r_tibia(2), r_tibia(3)-90];

    % display the femur cs on top of tibia cs (rotation become easier to  understand and compared)
    display_axis_v2(ax1(1), t_tibia, R_tmp * eye(3), 15, '', 'Tag', 'cs_bone_gt', 'Mode', 'shadow', 'HideVector', [1 0 0]);
    display_axis_v2(ax1(2), t_tibia, R_tmp * eye(3), 15, '', 'Tag', 'cs_bone_gt', 'Mode', 'shadow', 'HideVector', [0 0 1]);    
    display_axis_v2(ax1(3), t_tibia, R_tmp * eye(3), 15, '', 'Tag', 'cs_bone_gt', 'Mode', 'shadow', 'HideVector', [0 1 0]);

    % display the line from femur cs to tibia cs (translation become easier to understand and compared)
    line(ax1(1), [0 t_tibia(1)], [0 t_tibia(2)], [0 t_tibia(3)], 'LineStyle', '-', 'Color', 'm', 'LineWidth', 1, 'Tag', 'line_t_gt');
    line(ax1(2), [0 t_tibia(1)], [0 t_tibia(2)], [0 t_tibia(3)], 'LineStyle', '-', 'Color', 'm', 'LineWidth', 1, 'Tag', 'line_t_gt');
    line(ax1(3), [0 t_tibia(1)], [0 t_tibia(2)], [0 t_tibia(3)], 'LineStyle', '-', 'Color', 'm', 'LineWidth', 1, 'Tag', 'line_t_gt');


    % display the translation text
    text( ax1(1), t_tibia(1), t_tibia(2), t_tibia(3), sprintf('  (%.2f, %.2f)', t_tibia(2), t_tibia(3)), ...
          'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'Tag', 'text_t_gt');
    text( ax1(2), t_tibia(1), t_tibia(2), t_tibia(3), sprintf('(%.2f, %.2f)  ', t_tibia(1), t_tibia(3)), ...
          'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'Tag', 'text_t_gt');
    text( ax1(3), t_tibia(1), t_tibia(2), t_tibia(3), sprintf('(%.2f, %.2f)  ', t_tibia(1), t_tibia(2)), ...
          'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', 'Tag', 'text_t_gt');

    % display the rotation text
    scaler = 20;
    text( ax1(1), t_tibia(1), t_tibia(2)-15, t_tibia(3), sprintf('%.2f°', r_tibia(3)), ...
          'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'Tag', 'text_r_gt');
    text( ax1(2), t_tibia(1)+15, t_tibia(2), t_tibia(3), sprintf('%.2f°', r_tibia(2)), ...
          'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'Tag', 'text_r_gt');
    text( ax1(3), t_tibia(1), t_tibia(2)-15, t_tibia(3), sprintf('%.2f°', r_tibia(1)), ...
          'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'Tag', 'text_r_gt');

    % Test lagi -----------------------------------------------------------

    % loop for all dofs, let's display them one by one
    for idx_dof=1:6
        xline(ax2(idx_dof), idx_t, 'k', 'Tag', 'current_time');
        plot(ax2(idx_dof), idx_t, kneeJoint6DOFs_estSmooth(idx_t, idx_dof), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 10, 'Tag', 'current_est');
        plot(ax2(idx_dof), idx_t, kneeJoint6DOFs_gt(idx_t, idx_dof), 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 10, 'Tag', 'current_gt');


        % set the legend, but i only want the last two, that is my
        % estimation and ground truth line (they are stored in a way that 
        % the newest plot will be in the first element)
        mylines = findobj(ax2(idx_dof), 'Type','Line');
        legend(mylines(end-1: end), {'Ground truth', 'Estimation'}, 'FontSize', 12);


    end


    drawnow;

end