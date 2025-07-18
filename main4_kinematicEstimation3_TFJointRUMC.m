%% HELOOOOOO
% - Before running this script, you must already generated .mat file called
%   "all_TsReg_sX_mXX.mat". This .mat file contains all the T_bone_ref (the 
%   result of the regustration) for the whole timestamps. 
% - The purpose of this script is to process those Ts and convert it in the
%   knee joint coordinate system and quantify it against the ground truth.
% - The result of this script is a .mat file called "aall_TsGTEst_sX_mXX.mat"
%   which will be stored in a folder called "Tdata_sX_mXX_DATE-TIME". This 
%   .mat file contains transformations of the bones (femur and tibia) 
%   relative to the ref (global) both GT and est.

clc; clear; close all;

% [EDIT] directory to the project
path_root    = 'D:\Documents\BELANDA\PhD Thesis\Code\MATLAB\amode_navigation_experiment\experiment_b';

% [EDIT] directory to the trial
dir_trial    = "trial_0023_Session4_04";

% [EDIT] Change the data you are using accordingly. 
% ------ dir_depthdata is created by main1_processDepthData.m
% ------ dir_Tdata is created by main3_registrationWithTime.m
dir_depthdata = 'depthdata_s4_m04_20250708-172830';
dir_Tdata     = 'Tdata_s4_m04_20250710-074800';

% [EDIT] Change this to select the pin [femur, tibia], 1 -> PRO, 2-> DIS
idcs_pin = [2, 1];

% [EDIT] For displaying [amode3d, tibiaest, femur and tibiagt]
is_display.amode3d  = true;
is_display.tibiaest = true;
is_display.bonegt   = [true, true];

% [EDIT] for saving the resulting mat file
is_saveMat = false;


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

% declare the path which consists of processed depth data
path_outputdepth = fullfile(path_outputs, 'output_allest', dir_depthdata);

% Generate path to function directory
addpath(genpath(path_function));


%% INITIALIZE SOME DATA

% 1) Load and organize the struct for the CT data (allBone_CT)
run('extra_structCTdata.m');

% 2) Load the rigid body data that is stored in a .csv File
% -- Use the dir function to find all files with the .csv extension in the
%    directory. 
% -- This .csv file is the recording of the rigid bodies from Qualisys. 
% -- Here, we assume that there is always only one csv file the directory
% fileList = dir(fullfile(path_measurement, folders_measurement{folder_idx}, '*.csv'));
% 
% % Check if the .csv file was found
% if ~isempty(fileList)
%     % Assuming there's only one CSV file, get its full name
%     fullFileName = fullfile(path_measurement, folders_measurement{folder_idx}, fileList(1).name);
%     % Read the CSV file and get the rigid bodies data (Table object)
%     all_rigidbodies_table = readCSV_qualisysRigidBodies(fullFileName);
% else
%     % Display a message if no CSV file was found
%     disp('No CSV file found in the specified directory.');
%     return;
% end

% But because it is too long to load, i will use the shortcut (i generated
% the mat file already, check the snippet code for loading in
% main_2_process3Damode.m)
load('all_rigidbodies_table.mat');

% 3) Load the depth data (all_depthestmean_table and all_deptheststd_table)
tmp_str = split(dir_trial, '_');
sess_str = tmp_str{3};
meas_str = tmp_str{4};
mat_filename = sprintf('all_amode3d_s%s_m%s.mat', sess_str(end), meas_str);
mat_fullpath = fullfile(path_outputdepth, mat_filename);
load(mat_fullpath);

% 4) Load the T data from registration
mat_filename = sprintf('all_TsReg_s%s_m%s.mat', sess_str(end), meas_str);
mat_fullpath = fullfile(path_outputdepth, dir_Tdata, mat_filename);
load(mat_fullpath);

% 3) There is this thing called REF. 
% -- It is a calibration tool for B-mode ultrasound. We were using this for 
%    bone surface reconstruction as well.
% -- They were all represented in respect of this REF. So naturally we are 
%    using this REF as our reference for everything. 
% -- Well, it is actually not really necessary to do that, we can just use 
%    Qulisys coordinate frame, but because in the data processing in the 
%    earlier experiment we are using REF coordinate frame, so let's just 
%    use it again.
ref_name = 'B_N_REF';


%% INITIALIZE FIGURE OBJECTS AND EVERYTHING RELATED TO PLOTS

fig1 = figure('Name', 'Kinematic Estimation');
ax1 = axes(fig1);
axis(ax1, 'equal');
hold(ax1, 'on');
grid(ax1, 'on');
xlabel(ax1, 'X');
ylabel(ax1, 'Y');
zlabel(ax1, 'Z');
view(ax1, 30,15);


%% MAIN PROGRAM

% Initialize the first position of the bone
boneCTpoints_init = T_prereg_boneCT * cartesian2homogeneous(allBone_CT(2).stl.Points);

% Initialized the valid timestamps
timestamp_idcs_valid = all_TsReg_table.Timestamp_idx;
timestamp_ms_valid   = all_TsReg_table.Timestamp_ms;
n_timestamp_valid    = length(timestamp_idcs_valid);

% Allocate memory for storing all femur and tibia (est and gt)
% transformation. Later we will put them into one single table
Ts_femurGT_ref      = cell(n_timestamp_valid,1);
Ts_tibiaGT_ref      = cell(n_timestamp_valid,1);
Ts_tibiaEst_ref     = cell(n_timestamp_valid,1);
kneeJoint6DOFs_gt   = cell(n_timestamp_valid,1);
kneeJoint6DOFs_est  = cell(n_timestamp_valid,1);
errors_kneeJoint6DOF = cell(n_timestamp_valid, 1);

% Post process the Ts_boneEst_ref we have.
Ts_boneEst_ref_3dmat = cat(3, all_TsReg_table.Ts_boneEst_ref{:});
Ts_boneEstSmooth_ref_3dmat = smoothTransformations(Ts_boneEst_ref_3dmat, 'method', 'sgolay', 'window', 50);
Ts_boneEstSmooth_ref = squeeze( num2cell(Ts_boneEstSmooth_ref_3dmat, [1 2]) );
all_TsReg_table.Ts_boneEstSmooth_ref = Ts_boneEstSmooth_ref;

% initialize waitbar, so the user not get bored
h = waitbar(0, 'Please wait...');

% loop for all timeframe
for idx_t_3damode = 1:n_timestamp_valid

    % show the waitbar progress
    waitbar(idx_t_3damode/n_timestamp_valid, h, sprintf('Timestamp processed: %d/%d', idx_t_3damode, n_timestamp_valid));

    delete(findobj(ax1, 'Tag', '3d_amode'));
    delete(findobj(ax1, 'Tag', '3d_bone_est'));
    delete(findobj(ax1, 'Tag', '3d_bone_gt'));
    delete(findobj(ax1, 'Tag', 'cs_bone_est'));
    delete(findobj(ax1, 'Tag', 'cs_bone_gt'));

    % 1) AMODE 3D PART ----------------------------------------------------
    
    % 1.1) Get the current amode 3d points
    current_boneQUSpoints = all_amode3d_table{idx_t_3damode, 'Points'};
    current_boneQUSpoints = current_boneQUSpoints{1};

    % 1.2) Display the current amode 3d points
    if(is_display.amode3d)
        scatter3( ax1, ...
                  current_boneQUSpoints(1,:), ...
                  current_boneQUSpoints(2,:),  ...
                  current_boneQUSpoints(3,:), ...
                  50, 'red', 'filled', ...
                  'Tag', '3d_amode');
    end


    % 2) BONE CS ESTIMATION PART ------------------------------------------

    % 2.1) Transform the bone to it current 3d pose
    T_icpnow_icpinit = all_TsReg_table.Ts_icpnow_icpinit{idx_t_3damode};
    current_boneCTpoints_est = T_icpnow_icpinit * boneCTpoints_init;
    current_boneCTtri_est    = triangulation(allBone_CT(2).stl.ConnectivityList, current_boneCTpoints_est(1:3,:)');
    
    % 2.2) Get and store the current coordinate system
    T_boneEst_ref = all_TsReg_table.Ts_boneEstSmooth_ref{idx_t_3damode};
    Ts_tibiaEst_ref{idx_t_3damode} = T_boneEst_ref;

    % 2.4) Show the bone and the display axes
    if(is_display.tibiaest)
        trisurf(current_boneCTtri_est, 'FaceColor', '#2ecc71', 'FaceAlpha', 0.15, 'EdgeColor', '#2ecc71', 'EdgeAlpha', 0.1, 'Parent', ax1, 'Tag', '3d_bone_est')
        display_axis(ax1, T_boneEst_ref(1:3, 4), T_boneEst_ref(1:3, 1:3), 30, 'T_boneEst_ref', 'Tag', 'cs_bone_est');
    end

    % 3) GROUND TRUTH CS GT PART ------------------------------------------

    % Loop for two bones. The code to display is practically identical.
    for idx_bone = 1:2

        % 3.0) Select which pin should we use
        idx_pin = idcs_pin(idx_bone);

        % 3.1) Get the relevant transformation of the bone from CT scan
        T_bone_ct = allBone_CT(idx_bone).T;
        bone_name = allBone_CT(idx_bone).name;
        T_pin_ct  = allBone_CT(idx_bone).pin(idx_pin).T;
        pin_name  = allBone_CT(idx_bone).pin(idx_pin).name;
    
        % 3.2) Get the current T_pin_Q
        idx_t_rigidbody = all_amode3d_table{idx_t_3damode, 'Timestamp_idx'};
        T_pin_Q = all_rigidbodies_table.(pin_name)(idx_t_rigidbody);
        T_pin_Q = T_pin_Q{1}.T;
        T_ref_Q = all_rigidbodies_table.(ref_name)(idx_t_rigidbody);
        T_ref_Q = T_ref_Q{1}.T;
    
        % 3.3) Calculate transformation of pin in ref 
        T_pin_ref = inv(T_ref_Q) * T_pin_Q;
    
        % 3.4) Calculate transformation of boneGT in ref
        T_bone_pin    = inv(T_pin_ct) * T_bone_ct;
        T_boneGT_ref  = T_pin_ref * T_bone_pin;
    
        % 3.5) Transform the point cloud
        boneCTpoints_original   = cartesian2homogeneous(allBone_CT(idx_bone).stl.Points);
        current_boneCTpoints_GT = T_boneGT_ref * inv(T_bone_ct) * boneCTpoints_original;
        current_boneCTtri_GT    = triangulation(allBone_CT(idx_bone).stl.ConnectivityList, current_boneCTpoints_GT(1:3,:)');

        % 3.6) Store the femur and tibia bone CS gt in ref (global)
        if(idx_bone==1)
            Ts_femurGT_ref{idx_t_3damode} = T_boneGT_ref;
        else
            Ts_tibiaGT_ref{idx_t_3damode} = T_boneGT_ref;
        end

        % 3.7) Display the Bone
        if(is_display.bonegt(idx_bone))
            trisurf(current_boneCTtri_GT, 'FaceColor', '#bdc3c7', 'FaceAlpha', 0.15, 'EdgeColor', 'k', 'EdgeAlpha', 0.1, 'Parent', ax1, 'Tag', '3d_bone_gt');
            display_axis(ax1, T_boneGT_ref(1:3, 4), T_boneGT_ref(1:3, 1:3), 30, 'T_boneGT_ref', 'Tag', 'cs_bone_gt');
        end

    % end bone loop    
    end


    % 3) CALCULATE JOINT KINEMATICS (EST AND GT) --------------------------

    % 1) Calculate the knee joint 6Dof values (gt and est)
    % % % The function below is mine, i am not sure that it is correct or not
    [T_est, r_est, t_est] = generateKneeJointCS_v2(Ts_femurGT_ref{idx_t_3damode}, T_boneEst_ref, 'r');
    [T_gt, r_gt, t_gt]    = generateKneeJointCS_v2(Ts_femurGT_ref{idx_t_3damode}, Ts_tibiaGT_ref{idx_t_3damode}, 'r');
    r_est = rad2deg(r_est);
    r_gt = rad2deg(r_gt);

    % % The function below, i got from Miriam (ORL, RadboudUMC). They usually
    % % construct the R matrix by stacking the basis vector of the ACS
    % % row-wise (which is actually not standard, it should be column-wise).
    % % Consequently the function below also requires the similar structure. 
    % % Since i use the column-wise R, before i feed to this function i need
    % % to transpose it first.
    % T_femurGT_ref = Ts_femurGT_ref{idx_t_3damode};
    % r_est = findang_groodsuntay_style(T_femurGT_ref(1:3, 1:3)', T_boneEst_ref(1:3, 1:3)', 'right');
    % % This one, i have no clue why this operation, but i just following
    % % them. Note. They use row translation vector, i use column, so i need
    % % to transpose it first; See the R in the end of the operation? here i 
    % % didn't put a transpose, because the original one put a transpose on it.
    % t_est = (T_boneEst_ref(1:3, 4)' - T_femurGT_ref(1:3, 4)') * T_femurGT_ref(1:3, 1:3);
    % 
    % % Below is basically the same as above, just for GT
    % T_tibiaGT_ref = Ts_tibiaGT_ref{idx_t_3damode};
    % r_gt = findang_groodsuntay_style(T_femurGT_ref(1:3, 1:3)', T_tibiaGT_ref(1:3, 1:3)', 'right');
    % t_gt = (T_tibiaGT_ref(1:3, 4)' - T_femurGT_ref(1:3, 4)') * T_femurGT_ref(1:3, 1:3);


    % % My implementation of what Miriam did, which is basically just Tibia
    % % ACS relative to Femur ACS. Below is for estimation
    % T_femurGT_ref      = Ts_femurGT_ref{idx_t_3damode};
    % T_tibiaEst_femurGT = inv(T_femurGT_ref) * T_boneEst_ref;
    % r_est = rad2deg(rotm2eul(T_tibiaEst_femurGT(1:3, 1:3), 'ZYX'));
    % t_est = T_tibiaEst_femurGT(1:3, 4)';
    % % Below is for ground truth
    % T_tibiaGT_ref      = Ts_tibiaGT_ref{idx_t_3damode};
    % T_tibiaGT_femurGT = inv(T_femurGT_ref) * T_tibiaGT_ref;
    % r_gt = rad2deg(rotm2eul(T_tibiaGT_femurGT(1:3, 1:3), 'ZYX'));
    % t_gt = T_tibiaGT_femurGT(1:3, 4)';

    %  Store the GT
    kneeJoint6DOFs_est{idx_t_3damode}  = [t_est, r_est];
    kneeJoint6DOFs_gt{idx_t_3damode}   = [t_gt, r_gt];
    errors_kneeJoint6DOF{idx_t_3damode} = abs(kneeJoint6DOFs_est{idx_t_3damode} - kneeJoint6DOFs_gt{idx_t_3damode});

    % draw the plot
    drawnow;

% end timestamp loop
end

% Close the waitbar
close(h);

% Construct the table
all_kneeJoint6DOFs_table = table( timestamp_idcs_valid, timestamp_ms_valid, ...
                                  Ts_femurGT_ref, Ts_tibiaGT_ref, Ts_tibiaEst_ref, ...
                                  kneeJoint6DOFs_gt, kneeJoint6DOFs_est, errors_kneeJoint6DOF, ...
                                  'VariableNames', ...
                                  {'Timestamp_idx', 'Timestamp_ms', ...
                                  'Ts_femurGT_ref', 'Ts_tibiaGT_ref', 'Ts_tibiaEst_ref', ...
                                  'kneeJoint6DOFs_gt', 'kneeJoint6DOFs_est', 'error_kneeJoint6DOF'});

% Save the data
if (is_saveMat)
    % save
    mat_filename = sprintf('all_kneeJoint6DOFs_s%s_m%s.mat', sess_str(end), meas_str);
    mat_fullpath = fullfile(path_outputdepth, dir_Tdata, mat_filename);
    save(mat_fullpath, 'all_kneeJoint6DOFs_table');
end