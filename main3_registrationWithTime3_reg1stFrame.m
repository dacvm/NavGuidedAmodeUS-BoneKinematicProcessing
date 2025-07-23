%% HELOOOO
% - Before running this script, you must already generated .mat file called
%   "all_amode3d_sX_mXX.mat". This .mat file contains well-structured amode
%   3d point cloud complete with its label (holder group) for all timestamp.
% - The purpose of this script is to register the known bone model to those
%   amode 3d point cloud.
% - The result of this script is a .mat file called "all_TsReg_sX_mXX.mat"
%   which will be stored in a folder called "Tdata_sX_mXX_DATE-TIME". This
%   mat file contains all the T_bone_ref (the result of the regustration) 
%   for the whole timestamps. 

clc; clear; close all;

% [EDIT] directory to the project
path_root    = 'D:\DennisChristie\NavGuidedAmodeUS-BoneKinematicProcessing';

% [EDIT] directory to the trial
% dir_trial    = "trial_0023_Session4_04";
dir_trial    = "trial_0011_Session3_02";

% [EDIT]
% dir_depthdata = 'depthdata_s4_m04_20250708-172830';
% dir_depthdata = 'depthdata_s3_m02_20250722-114731';
dir_depthdata = 'depthdata_s3_m02_20250722-174503';

% [EDIT] Select bone and pin
idx_bone = 2;
idx_pin  = 1;

% [EDIT]
is_displayRegProcess = false;
is_saveMat = true;

% [EDIT]
% is_usenavigationdata = true;

% [EDIT] Everything related to PICP
% path
path_picp    = "D:\DennisChristie\SwarmPerturbation-ICP";
% parameters
params_picp.name               = 'tibia';
params_picp.max_iters          = 10;
params_picp.rmse_threshold     = 0.001;
params_picp.init_perturb_rot   = 1.0;
params_picp.init_perturb_trans = 1.0;
params_picp.decay_rate         = 0.01;
params_picp.n_candidate        = 16;

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
addpath(path_picp);


% 6) Create folder for output if user want to save
if (is_saveMat)
    % get the session name
    tmp_str = split(dir_trial, '_');
    sess_str = tmp_str{3};
    meas_str = tmp_str{4};

    % get the current time
    foldername_timestamp  = datestr(now, 'yyyymmdd-HHMMSS');
    % create a unique name
    foldername_withprefix = sprintf('Tdata_s%s_m%s_%s', sess_str(end), meas_str, foldername_timestamp);
    % Full path and creation
    mkdir(fullfile(path_outputdepth, foldername_withprefix));
end


%% INITIALIZE SOME RIGID BODY METADATA

% 1) Create base rotation to trasnform qualisys base vector to matlab
% -- Qualisys has y direction as up, MATLAB has z direction as up
R_tmp = eul2rotm([0 0 pi/2], "ZYX");
t_tmp = [0 0 0]';
baseRotation_Qualisys2Matlab = [R_tmp, t_tmp; 0 0 0 1];

% 2.a) Load and organize the struct for the CT data (allBone_CT)
run('extra_structCTdata.m');
run('extra_structPrereg.m')

% 2.b) Get the relevant point cloud of the bone from CT scan
boneCTstl_original     = allBone_CT(idx_bone).stl;
boneCTpc_original      = pointCloud(boneCTstl_original.Points);
boneCTnormals_original = STLVertexNormals(boneCTstl_original.ConnectivityList, boneCTstl_original.Points)';
boneCTpoints_original  = [boneCTpc_original.Location'; ones(1, length(boneCTpc_original.Location))];

% 2.c) Get the relevant transformation of the bone from CT scan
T_bone_ct = allBone_CT(idx_bone).T;
bone_name = allBone_CT(idx_bone).name;
T_pin_ct  = allBone_CT(idx_bone).pin(idx_pin).T;
pin_name  = allBone_CT(idx_bone).pin(idx_pin).name;

% 3) There is this thing called REF. 
% -- It is a calibration tool for B-mode ultrasound. We were using this for 
%    bone surface reconstruction as well.
% -- They were all represented in respect of this REF. So naturally we are 
%    using this REF as our reference for everything. 
% -- Well, it is actually not really necessary to do that, we can just use 
%    Qulisys coordinate frame, but because in the data processing in the 
%    earlier experiment we are using REF coordinate frame, so let's just 
%    use it again.
ref_name   = 'B_N_REF';

%% INITIALIZE FIGURE OBJECTS AND EVERYTHING RELATED TO PLOTS

if(is_displayRegProcess)
    fig1 = figure('Name', 'Registration Process');
    ax1  = axes(fig1);
    axis(ax1, 'equal');
    hold(ax1, 'on');
    grid(ax1, 'on');
    xlabel(ax1, 'X');
    ylabel(ax1, 'Y');
    zlabel(ax1, 'Z');
    view(ax1, 30,15);
end


%% INITIAL-STAGE: PRE-REGISTRATION

% 1) Load the rigid body data that is stored in a .csv File
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
% load('all_rigidbodies_table_s4_m04_d01.mat');
load('all_rigidbodies_table_s3_m02_d01.mat');

% 2) Load the depth data (all_depthestmean_table and all_deptheststd_table)
tmp_str = split(dir_trial, '_');
sess_str = tmp_str{3};
meas_str = tmp_str{4};
mat_filename = sprintf('all_amode3d_s%s_m%s.mat', sess_str(end), meas_str);
load(fullfile(path_outputdepth, mat_filename));

% 3) Get the initial data for preregistration purposes
boneQUSpoints_ref = all_amode3d_table{1, 'Points'};
boneQUSpoints_ref = boneQUSpoints_ref{1};
boneQUSamps   = all_amode3d_table{1, 'Ampitude'};
boneQUSamps   = boneQUSamps{1};
boneQUSlabels = all_amode3d_table{1, 'Label'};
boneQUSlabels = boneQUSlabels{1};

% 4) Select the bone for prereg area
bonePrereg = allBone_preReg(idx_bone);

% 5) Perform preregistration. This is basically similar script as in the 
% previous experiment (experiment_a)
[T_prereg_boneCT, boneCTpoints_prereg] = preRegistration( boneCTpoints_original, boneQUSpoints_ref, boneQUSlabels, bonePrereg, false );

%%

% This matrix will store transformation of current timestamp registration 
% relative to the previous timestamp registration
T_icpnow_icpprev = eye(4);
% This matrix will store accumulative transformation
T_icpnow_icpinit = eye(4);

% Renaming
current_boneCTpoints_est = boneCTpoints_prereg;
init_boneCTpoints_est    = boneCTpoints_prereg;

% Initialized the valid timestamps
timestamp_idcs_valid = all_amode3d_table.Timestamp_idx;
timestamp_ms_valid   = all_amode3d_table.Timestamp_ms;
n_timestamp_valid    = length(timestamp_idcs_valid);

% Dummy for debugging purposes only
n_timestamp_valid = 1310;
timestamp_idcs_valid = timestamp_idcs_valid(1:n_timestamp_valid);
timestamp_ms_valid   = timestamp_ms_valid(1:n_timestamp_valid);

% allocate memory for storing all registration transformation. Later we
% will put them into one single table
Ts_icpnow_icpprev   = cell(n_timestamp_valid,1);
Ts_icpnow_icpinit   = cell(n_timestamp_valid,1);
Ts_icpnowKF_icpinit = cell(n_timestamp_valid,1);
Ts_boneEst_ref      = cell(n_timestamp_valid,1);
Ts_boneEstKF_ref    = cell(n_timestamp_valid,1);
Ts_boneGT_ref       = cell(n_timestamp_valid,1);
errors_T                   = cell(n_timestamp_valid, 1);
errors_rmse_us2regbone     = zeros(n_timestamp_valid, 1);
errors_rmse_regbone2gtbone = zeros(n_timestamp_valid, 1);

% initialize waitbar, so the user not get bored
h = waitbar(0, 'Please wait...');

% test
ukf  = [];

% loop for all timeframe
for idx_t_3damode = 1:n_timestamp_valid

    % show the waitbar progress
    waitbar(idx_t_3damode/n_timestamp_valid, h, sprintf('Timestamp processed: %d/%d', idx_t_3damode, n_timestamp_valid));

    % 0.a) Delete objects in the figure
    if(is_displayRegProcess)
        delete(findobj(ax1, 'Tag', '3d_amode'));
        delete(findobj(ax1, 'Tag', '3d_bone_est'));
        delete(findobj(ax1, 'Tag', '3d_bone_gt'));
        delete(findobj(ax1, 'Tag', 'cs_bone_est'));
        delete(findobj(ax1, 'Tag', 'cs_bone_gt'));
    end
    
    % 0.b) Start the registration process
    fprintf('Registeration t=%d ', idx_t_3damode); tic;


    % 1) Registration Estimation Part -------------------------------------

    % 1.1) Get the initial data for preregistration purposes
    current_boneQUSpoints = all_amode3d_table{idx_t_3damode, 'Points'};
    current_boneQUSpoints = current_boneQUSpoints{1};

    % 1.2) Register with icp(moving, fixed)
    [T_currentBoneQUS_currentBoneCT, ~] = icp_with_perturbation_v2( current_boneCTpoints_est, current_boneQUSpoints, ...
                                                                    params_picp.max_iters, ...
                                                                    params_picp.rmse_threshold, ...
                                                                    params_picp.init_perturb_rot, ...
                                                                    params_picp.init_perturb_trans, ...
                                                                    params_picp.decay_rate, ...
                                                                    params_picp.n_candidate, ...
                                                                    false, false, false);

    % - I rename it, so Dennis in the future understand the context:
    % - Technically, this is a transformation that transform the current 
    %   timestamp of boneCT to boneQUS (or boneQUS in respect to BoneCT).
    % - But in the context of consecutive timestamps, this technical name
    %   of the T seems so vague, so, instead, i use this naming system.
    T_icpnow_icpprev = T_currentBoneQUS_currentBoneCT;

    % 1.3a) Accumulate the T
    T_icpnow_icpinit = T_icpnow_icpprev * T_icpnow_icpinit;

    % Kalman Filter?
    if(idx_t_3damode>3)
        Ts_icpnow_icpinit_cell    = Ts_icpnow_icpinit(idx_t_3damode-3:idx_t_3damode-1);
        Ts_icpnow_icpinit_3dmat   = cat(3, Ts_icpnow_icpinit_cell{:});

        kalmanposese3.Ts_est      = Ts_icpnow_icpinit_3dmat;
        kalmanposese3.t           = timestamp_ms_valid(idx_t_3damode-3:idx_t_3damode-1);
        kalmanposese3.T_meas_i    = T_icpnow_icpinit;
        kalmanposese3.t_i         = timestamp_ms_valid(idx_t_3damode);

        if(idx_t_3damode == 895)
            a=1;
        end

        [T_icpnow_icpinit, ukf] = kalmanPoseSE3_v2(kalmanposese3.Ts_est, ...
                                                   kalmanposese3.t, ...
                                                   kalmanposese3.T_meas_i, ...
                                                   kalmanposese3.t_i, ...
                                                   ukf);
    end

    % Transform the bone
    current_boneCTpoints_est = T_icpnow_icpinit * init_boneCTpoints_est;

    % 1.3b) Calculate the T in ref. 
    % ----> This will be usefull for calculating the T error later
    T_boneEst_ref    = T_icpnow_icpinit * T_prereg_boneCT * T_bone_ct;

    % 1.3c) Construct the triangulation object for the registered bone. 
    % ----> This will be usefull for calculating the RMSE later
    current_boneCTtri_est    = triangulation(boneCTstl_original.ConnectivityList, current_boneCTpoints_est(1:3,:)');

    % 1.4) Store the transformation
    % Ts_icpnow_icpprev{idx_t_3damode} = T_icpnow_icpprev;
    Ts_icpnow_icpinit{idx_t_3damode}   = T_icpnow_icpinit;
    Ts_boneEst_ref{idx_t_3damode}      = T_boneEst_ref;


    % 2) Ground Truth Part ------------------------------------------------

    % 2.1) Get the current T_pin_Q
    idx_t_rigidbody = all_amode3d_table{idx_t_3damode, 'Timestamp_idx'};
    T_pin_Q = all_rigidbodies_table.(pin_name)(idx_t_rigidbody);
    T_pin_Q = T_pin_Q{1}.T;
    T_ref_Q = all_rigidbodies_table.(ref_name)(idx_t_rigidbody);
    T_ref_Q = T_ref_Q{1}.T;

    % 2.2) Calculate transformation of pin in ref 
    T_pin_ref = inv(T_ref_Q) * T_pin_Q;

    % 3.3) Calculate transformation of boneGT in ref
    T_bone_pin    = inv(T_pin_ct) * T_bone_ct;
    T_boneGT_ref  = T_pin_ref * T_bone_pin;

    % 2.4) Transform the point cloud
    current_boneCTpoints_GT = T_boneGT_ref * inv(T_bone_ct) * boneCTpoints_original;
    current_boneCTtri_GT    = triangulation(boneCTstl_original.ConnectivityList, current_boneCTpoints_GT(1:3,:)');

    % 2.5) Store the transformation
    Ts_boneGT_ref{idx_t_3damode} = T_boneGT_ref;



    % 3) Error Calculation ------------------------------------------------

    % 3.1) Calculate the transformation error
    [t_err, r_err] =  calculateRBerror(T_boneGT_ref, T_boneEst_ref);
    fprintf('| T_est_err = [%.2f, %.2f, %.2f, %.2f, %.2f, %.2f]\t', [t_err, r_err]);


    % 3.2) Calculate the mean distance of amode 3d points to the registered bone
    meanError1_pt2plane = calculateRMSE(current_boneQUSpoints(1:3, :)', current_boneCTtri_est, 'pt2plane');

    % 3.3) Calculate the mean distance of registered bone to gt bone
    meanError2_pt2plane = calculateRMSE(current_boneCTtri_est,  current_boneCTtri_GT, 'pt2plane');

    % 3.4) Store the errors
    errors_T{idx_t_3damode} = [t_err, r_err];
    errors_rmse_us2regbone(idx_t_3damode)     = meanError1_pt2plane;
    errors_rmse_regbone2gtbone(idx_t_3damode) = meanError2_pt2plane;


    % 4) Display and Others -----------------------------------------------

    if(is_displayRegProcess)
        % Display the amode 3d
        scatter3( ax1, ...
                  current_boneQUSpoints(1,:), ...
                  current_boneQUSpoints(2,:),  ...
                  current_boneQUSpoints(3,:), ...
                  50, 'red', 'filled', ...
                  'Tag', '3d_amode');    

        % Display the registered bone
        display_axis(ax1, T_boneEst_ref(1:3, 4), T_boneEst_ref(1:3, 1:3), 30, 'T_boneEst_ref', 'Tag', 'cs_bone_est');
        trisurf(current_boneCTtri_est, 'FaceColor', '#2ecc71', 'FaceAlpha', 0.15, 'EdgeColor', '#2ecc71', 'EdgeAlpha', 0.1, 'Parent', ax1, 'Tag', '3d_bone_est')

        % Display the GT bone
        display_axis(ax1, T_boneGT_ref(1:3, 4), T_boneGT_ref(1:3, 1:3), 30, 'T_boneGT_ref', 'Tag', 'cs_bone_gt');
        trisurf(current_boneCTtri_GT, 'FaceColor', '#bdc3c7', 'FaceAlpha', 0.15, 'EdgeColor', 'k', 'EdgeAlpha', 0.1, 'Parent', ax1, 'Tag', '3d_bone_est');

        % draw to the plot
        drawnow;
    end


    % Display duration
    time = toc;
    fprintf('(%.2fs)\n', time);

end

% Close the waitbar
close(h);

%%


% Ts_boneEst_ref_3dmat       = cat(3, Ts_boneEst_ref{:});
% ts_boneEst_ref             = squeeze(Ts_boneEst_ref_3dmat(1:3, 4, :));
% Ts_boneEstSmooth_ref_3dmat = smoothTransformations(Ts_boneEst_ref_3dmat, 'method', 'sgolay', 'window', 30);
% ts_boneEstSmooth_ref       = squeeze(Ts_boneEstSmooth_ref_3dmat(1:3, 4, :));
% 
% Ts_boneGT_ref_3dmat    = cat(3, Ts_boneGT_ref{:});
% ts_boneGT_ref = squeeze(Ts_boneGT_ref_3dmat(1:3, 4, :));
% 
% fig2 = figure;
% ax2 = axes(fig2);
% plot3(ax2, ts_boneGT_ref(1,:), ts_boneGT_ref(2,:), ts_boneGT_ref(3,:), '-ob'); grid on; axis equal; hold on;
% plot3(ax2, ts_boneEstSmooth_ref(1,:), ts_boneEstSmooth_ref(2,:), ts_boneEstSmooth_ref(3,:), '-or');
% 
% errors_T_mat  = cat(1, errors_T{:});
% errors_T_mean = mean(abs(errors_T_mat), 1);
% fprintf('Mean Absolute Error T = %.2f, %.2f, %.2f, %.2f, %.2f, %.2f\n', errors_T_mean);

% Post process the Ts_boneEst_ref we have.
Ts_boneEst_ref_3dmat = cat(3, Ts_boneEst_ref{:});
Ts_boneEstSmooth_ref_3dmat = smoothTransformations(Ts_boneEst_ref_3dmat, 'method', 'sgolay', 'window', 30);
Ts_boneEstSmooth_ref = squeeze( num2cell(Ts_boneEstSmooth_ref_3dmat, [1 2]) );

% Construct the table
all_TsReg_table = table( timestamp_idcs_valid, timestamp_ms_valid, ...
                         Ts_icpnow_icpprev, Ts_icpnow_icpinit, Ts_boneEst_ref, Ts_boneEstSmooth_ref, Ts_boneGT_ref, ...
                         errors_T, errors_rmse_us2regbone, errors_rmse_regbone2gtbone, ...
                         'VariableNames', ...
                         {'Timestamp_idx', 'Timestamp_ms', ...
                         'Ts_icpnow_icpprev', 'Ts_icpnow_icpinit', 'Ts_boneEst_ref', 'Ts_boneEstSmooth_ref', 'Ts_boneGT_ref', ...
                         'errors_T', 'errors_rmse_us2regbone', 'errors_rmse_regbone2gtbone'});


% Save the data
if (is_saveMat)
    % save
    mat_filename = sprintf('all_TsReg_s%s_m%s.mat', sess_str(end), meas_str);
    mat_fullpath = fullfile(path_outputdepth, foldername_withprefix, mat_filename);
    save(mat_fullpath, 'T_prereg_boneCT', 'all_TsReg_table');
end