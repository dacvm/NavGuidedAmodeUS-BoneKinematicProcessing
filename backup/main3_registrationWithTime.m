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
path_root    = 'D:\Documents\BELANDA\PhD Thesis\Code\MATLAB\amode_navigation_experiment\experiment_b';

% [EDIT] directory to the trial
dir_trial    = "trial_0023_Session4_04";

% [EDIT]
dir_depthdata = 'depthdata_s4_m04_20250708-172830';

% [EDIT] Select bone and pin
idx_bone = 2;
idx_pin  = 1;

% [EDIT]
is_displayRegProcess = false;
is_saveMat = true;

% [EDIT]
is_usenavigationdata = true;

% [EDIT] Everything related to PICP
% path
path_picp    = "D:\Documents\MATLAB\icp_with_pertubation";
% parameter
params_picp.name               = 'tibia';
params_picp.max_iters          = 50;
params_picp.rmse_threshold     = 0.001;
params_picp.init_perturb_rot   = 1.0;
params_picp.init_perturb_trans = 1.0;
params_picp.decay_rate         = 0.01;
params_picp.n_candidate        = 64;

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
load('all_rigidbodies_table.mat');

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
[T_prereg_bonect, ~] = preRegistration( boneCTpoints_original, boneQUSpoints_ref, boneQUSlabels, bonePrereg, false );

% 6) initialize T
boneCTpoints_prereg = T_prereg_bonect * boneCTpoints_original;



%% MAIN PROGRAM

% Initialize some data
current_boneCTpoints_est = boneCTpoints_prereg;
current_T = eye(4);

% Initialized the valid timestamps
timestamp_idcs_valid = all_amode3d_table.Timestamp_idx;
% n_timestamp_valid = length(timestamp_idcs_valid);

% Dummy for debugging purposes only
n_timestamp_valid = 1000;
timestamp_idcs_valid = timestamp_idcs_valid(1:n_timestamp_valid);

% allocate memory for storing all registration transformation. Later we
% will put them into one single table
% Ts_current_init = repmat(eye(4), 1, 1, n_timestamp_valid);
% Ts_current_prev = repmat(eye(4), 1, 1, n_timestamp_valid);
Ts_current_init   = cell(n_timestamp_valid,1);
Ts_current_prev   = cell(n_timestamp_valid,1);

% loop for all timeframe
for idx_t_3damode = 1:n_timestamp_valid

    % 0.a) Delete objects in the figure
    if(is_displayRegProcess)
        delete(findobj(ax1, 'Tag', 'amode_3d'));
        delete(findobj(ax1, 'Tag', 'bone_3d'));
    end
    
    % 0.b) Start the registration process
    fprintf('Registeration t=%d ', idx_t_3damode); tic;

    % 1) Get the initial data for preregistration purposes
    current_boneQUSpoints = all_amode3d_table{idx_t_3damode, 'Points'};
    current_boneQUSpoints = current_boneQUSpoints{1};
    current_boneQUSamps   = all_amode3d_table{idx_t_3damode, 'Ampitude'};
    current_boneQUSamps   = current_boneQUSamps{1};

    % 2) Register with icp(moving, fixed)
    [T_currentBoneQUS_currentBoneCT, ~] = icp_with_perturbation_v2( current_boneCTpoints_est, current_boneQUSpoints, ...
                                                                   params_picp.max_iters, ...
                                                                   params_picp.rmse_threshold, ...
                                                                   params_picp.init_perturb_rot, ...
                                                                   params_picp.init_perturb_trans, ...
                                                                   params_picp.decay_rate, ...
                                                                   params_picp.n_candidate, ...
                                                                   false, false, false);

    % 3) Update the 3d pose of the bone
    current_boneCTpoints_est = T_currentBoneQUS_currentBoneCT * current_boneCTpoints_est;
    current_T = T_currentBoneQUS_currentBoneCT * current_T;

    % 4) Store the transformation
    Ts_current_init{idx_t_3damode} = current_T;
    Ts_current_prev{idx_t_3damode} = T_currentBoneQUS_currentBoneCT;

    % 5) Display
    if(is_displayRegProcess)
        % 5.1) Display the amode 3d
        scatter3( ax1, ...
                  current_boneQUSpoints(1,:), ...
                  current_boneQUSpoints(2,:),  ...
                  current_boneQUSpoints(3,:), ...
                  50, 'red', 'filled', ...
                  'Tag', 'amode_3d');    

        % 5.2) Display the registered bone
        current_boneCTtri_est  = triangulation(boneCTstl_original.ConnectivityList, current_boneCTpoints_est(1:3,:)');
        trisurf(current_boneCTtri_est, 'FaceColor', '#2ecc71', 'FaceAlpha', 0.15, 'EdgeColor', '#2ecc71', 'EdgeAlpha', 0.1, 'Parent', ax1, 'Tag', 'bone_3d')

        drawnow;
    end

    % Display duration
    time = toc;
    fprintf('(%.2fs)\n', time);

end

% Construct the table
all_TsReg_table = table( timestamp_idcs_valid, Ts_current_prev, Ts_current_init, ...
                         'VariableNames', {'Timestamp_idx', 'Ts_current_prev', 'Ts_current_init'});


% Save the data
if (is_saveMat)
    % save
    mat_filename = sprintf('all_TsReg_s%s_m%s.mat', sess_str(end), meas_str);
    mat_fullpath = fullfile(path_outputdepth, foldername_withprefix, mat_filename);
    save(mat_fullpath, 'T_prereg_bonect', 'all_TsReg_table');
end






