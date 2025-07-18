%% HELOOOO
% Before running this script, you must already generated .mat file called
% "ust_depthtable.mat". This .mat file contains estimated bone depth for
% all ultrasound transducer for all timeframe. 

clc; clear; close all;

% [EDIT] directory to the project
path_root    = 'D:\Documents\BELANDA\PhD Thesis\Code\MATLAB\amode_navigation_experiment\experiment_b';

% [EDIT] directory to the trial
dir_trial    = "trial_0023_Session4_04";

% [EDIT] window configuration file
% csvfile_windowconfig = 'transducerconfig_v8a_window2024-12-20_14-37-59_edited2025-04-10_15-53-56.csv';

% [EDIT] holder configuration file
csvfile_holderconfig = 'transducerconfig_v8a.csv';

% [EDIT] Make sure you are using the correct depth data
dir_depthdata = 'depthdata_s4_m04_20250714-112513';

% [EDIT] Specify folder index
folder_idx = 1;

% [EDIT] Select bone and pin
idx_bone = 2;
idx_pin  = 1;

% [EDIT]
is_display = false;
is_saveMat = true;

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

% 1) Generate path to function directory
addpath(genpath(path_function));

% 2) Get all the folders of the intermediates
items_dir = dir(path_intermediate);
folders_intermediate = {items_dir([items_dir.isdir] & ~ismember({items_dir.name}, {'.', '..'})).name};

% 3) Get all the folders of the measurements
items_dir = dir(path_measurement);
folders_measurement = {items_dir([items_dir.isdir] & ~ismember({items_dir.name}, {'.', '..'})).name};
% the last folder for measurements is always empty (it was like that from
% the navigation system) we don't need that.
folders_measurement(end) = [];

% % 4) Load window configuration file
% csv_fullpath = fullfile(path_intermediate, csvfile_windowconfig);
% ust_windowconfig = table2struct(readtable(csv_fullpath));

% 5) Load holder configuration file
csv_fullpath = fullfile(path_data, 'configs', csvfile_holderconfig);
ust_holderconfig = table2struct(readtable(csv_fullpath));
% Find indices of structs where Group is not equal to 0, means that we are
% not using the ust
indices = arrayfun(@(x) x.Group ~= 0, ust_holderconfig);
% Create a new array of structs excluding those with Group equal to 0
ust_holderconfig = ust_holderconfig(indices);


%% INITIALIZE SOME RIGID BODY METADATA

% 1) Create base rotation to trasnform qualisys base vector to matlab
% -- Qualisys has y direction as up, MATLAB has z direction as up
R_tmp = eul2rotm([0 0 pi/2], "ZYX");
t_tmp = [0 0 0]';
baseRotation_Qualisys2Matlab = [R_tmp, t_tmp; 0 0 0 1];

% 2.a) Load and organize the struct for the CT data (allBone_CT)
run('extra_structCTdata.m');

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

if(is_display)
    fig1 = figure('Name', 'UST in 3d');
    ax1  = axes(fig1);
    axis(ax1, 'equal');
    hold(ax1, 'on');
    grid(ax1, 'on');
    xlabel(ax1, 'X');
    ylabel(ax1, 'Y');
    zlabel(ax1, 'Z');
    view(ax1, 30,15);
end


%% INITIALIZE SOME DATA

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

% 2) But because it is too long to load, i will use the shortcut (i generated
% the mat file already, check the snippet code for loading in
% main_2_process3Damode.m)
load('all_rigidbodies_table.mat');

% Load the depth data (all_depthestmean_table and all_deptheststd_table)
tmp_str = split(dir_trial, '_');
sess_str = tmp_str{3};
meas_str = tmp_str{4};
mat_filename = sprintf('all_depthest_s%s_m%s.mat', sess_str(end), meas_str);
load(fullfile(path_outputdepth, mat_filename));


%% MAIN PROGRAM

% Obtain valid timestamp
% -- Due to the nature of our depth detection algorithm (using bisplane), 
%    they extrapolate some parts of the depth in the timeline outside of
%    the detection data point (shoots up or down)
% -- We deleted that part by assigning bunch of zeros
% -- Each UST have different start and end of non-zeros, so here we want to
%    know those two indices.
tmp = table2array(all_depthestmean_table);
tmp = prod(tmp,2);
timestamp_idcs_valid = find(tmp~=0, 1, 'first'):find(tmp~=0, 1, 'last');
timestamp_ms_valid   = all_depthestmean_table.Timestamp_ms(timestamp_idcs_valid);

% Some initial variables
n_ust_selected  = size(all_depthestmean_table, 2);
n_timestamp     = size(timestamp_idcs_valid, 2);
ust_numbers_str = all_depthestmean_table.Properties.VariableNames(3:end);
ust_numbers     = str2double(ust_numbers_str);

% Variable to sore the resulting value
all_amode3d_matrix = zeros(n_timestamp, 3);

% Preallocate three N×1 cell arrays
cell_points = cell(length(timestamp_idcs_valid),1);
cell_amps   = cell(length(timestamp_idcs_valid),1);
cell_labels = cell(length(timestamp_idcs_valid),1);
idx_looptime = 1; 


% initialize waitbar, so the user not get bored
h = waitbar(0, 'Please wait...');

% loop only for valid timestamp
for idx_t = timestamp_idcs_valid
    %% LOOP FOR TIMESTAMPS

    % show the waitbar progress
    waitbar(idx_looptime/length(timestamp_idcs_valid), h, sprintf('Timestamp processed: %d/%d', idx_looptime, length(timestamp_idcs_valid)));

    % [0] -----------------------------------------------------------------
    % Delete every object in the plot before we go
    delete(findobj('Tag', 'cs_bonepin'));
    delete(findobj('Tag', 'cs_bonegt'));
    delete(findobj('Tag', 'cs_ust'));
    delete(findobj('Tag', 'bonesurface_ct'));
    delete(findobj('Tag', 'amode_3d'));

    
    % [1] -----------------------------------------------------------------
    % Get the current T_pin_Q
    T_pin_Q = all_rigidbodies_table.(pin_name)(idx_t);
    T_pin_Q = T_pin_Q{1}.T;
    T_ref_Q = all_rigidbodies_table.(ref_name)(idx_t);
    T_ref_Q = T_ref_Q{1}.T;


    % [2] -----------------------------------------------------------------
    % Calculate transformation of pin in ref 
    T_pin_ref = inv(T_ref_Q) * T_pin_Q;

    
    % [3] -----------------------------------------------------------------
    % Calculate transformation of boneGT in ref
    T_bone_pin    = inv(T_pin_ct) * T_bone_ct;
    T_boneGT_ref  = T_pin_ref * T_bone_pin;


    % [4] -----------------------------------------------------------------
    % Transform the point cloud
    T_ct_bone        = inv(T_bone_ct);
    boneCTpoints_GT  = T_boneGT_ref           * inv(T_bone_ct)      * boneCTpoints_original;
    boneCTnormals_GT = T_boneGT_ref(1:3, 1:3) * T_ct_bone(1:3, 1:3) * boneCTnormals_original;
    boneCTtri_GT     = triangulation(boneCTstl_original.ConnectivityList, boneCTpoints_GT(1:3,:)');

    % [5] -----------------------------------------------------------------
    % Display
    if(is_display)
        display_axis(ax1, T_pin_ref(1:3, 4), T_pin_ref(1:3, 1:3), 30, 'T_pin_ref', 'Tag', 'cs_bonepin');
        display_axis(ax1, T_boneGT_ref(1:3, 4), T_boneGT_ref(1:3, 1:3), 30, 'T_boneGT_ref', 'Tag', 'cs_bonegt');
        % trisurf(boneCTtri_GT, 'FaceColor', '#bdc3c7', 'FaceAlpha', 0.15, 'EdgeColor', 'k', 'EdgeAlpha', 0.1, 'Tag', 'bonesurface_ct');
    end

    % Preparing vars
    amode3d_points = zeros(4, length(ust_numbers));
    amode3d_amps   = zeros(length(ust_numbers), 1);
    amode3d_labels = cell(length(ust_numbers), 1);
    idx_loopust       = 1;

    % loop for all ust
    for idx_ust = ust_numbers
        %% LOOP FOR USTS
    
        % [1] -------------------------------------------------------------
        % Get the current ust and its group name
        current_uststruct    = ust_holderconfig(idx_ust);
        current_ustnumber    = ust_holderconfig(idx_ust).Number;
        current_ustgroupname = ust_holderconfig(idx_ust).GroupName;
        
        % [2] -------------------------------------------------------------
        % Transform the euler to rotation matrix
        % -- Here i used the equivalent lines as in the Qt C++ code. 
        % -- Why i do this, so that i know i am doing it exactly like i am 
        %    doing in the Qt C++ code (check VolumeAmodeVisualizer::updateTransformations)
        local_r_euler  = [current_uststruct.Local_Rx, current_uststruct.Local_Ry, current_uststruct.Local_Rz];
        local_R_matrix = rotx(local_r_euler(1)) * ...   % rotate about X by rx
                         roty(local_r_euler(2)) * ...   %   then about Y by ry
                         rotz(local_r_euler(3));        %   then about Z by rz
        % -- Get the translation vector
        local_t = [current_uststruct.Local_tx, current_uststruct.Local_ty, current_uststruct.Local_tz];
        % -- Create the T matrix
        T_probe_holder = [local_R_matrix, local_t'; 0 0 0 1];

        
        % [3] -------------------------------------------------------------
        % Get the current T_holder_Q for relevant idx_ust 
        current_RBstruct = all_rigidbodies_table.(current_ustgroupname)(idx_t);
        T_holder_Q       = current_RBstruct{1}.T;
        % -- Transform the Probe to REF coordinate frame
        T_probe_Q        = inv(T_ref_Q) * T_holder_Q * T_probe_holder;


        % [4] -------------------------------------------------------------
        % Get the current depth value for relevant idx_ust and idx_t 
        current_ustdepth         = all_depthestmean_table.(num2str(current_ustnumber))(idx_t);
        current_ustdepth3d_inQ   = [0 current_ustdepth 0 1]';
        current_ustdepth3d_inRef = T_probe_Q * baseRotation_Qualisys2Matlab * current_ustdepth3d_inQ;

        % Get the current amplitude value for relevant idx_ust and idx_t
        current_ustamp   = all_depthestamp_table.(num2str(current_ustnumber))(idx_t);


        % [5] -------------------------------------------------------------
        % Display
        if(is_display)
            axisname = sprintf('%d_%d', current_uststruct.Group, current_uststruct.Number);
            display_axis(ax1, T_probe_Q(1:3, 4), T_probe_Q(1:3, 1:3), 10, axisname, 'Tag', 'cs_ust');
            scatter3( ax1, ...
                      current_ustdepth3d_inRef(1), ...
                      current_ustdepth3d_inRef(2), ...
                      current_ustdepth3d_inRef(3), ...
                      50, 'red', 'filled', ...
                      'Tag', 'amode_3d');
        end


        % Collect necessary data from each ust
        amode3d_points(:,idx_loopust) = current_ustdepth3d_inRef;
        amode3d_amps(idx_loopust)     = current_ustamp;
        amode3d_labels{idx_loopust}   = current_ustgroupname;

        idx_loopust = idx_loopust+1;

    % end ust loop
    end
    
    % Store the depth mean and std to the big matrix
    cell_points{idx_looptime} = amode3d_points;
    cell_amps{idx_looptime}   = amode3d_amps;
    cell_labels{idx_looptime} = amode3d_labels;

    idx_looptime = idx_looptime+1;

% end timestamp loop
end

% Close the waitbar
close(h);

% Construct the table
all_amode3d_table = table( timestamp_idcs_valid', timestamp_ms_valid, cell_points, cell_amps, cell_labels, ...
                           'VariableNames', {'Timestamp_idx', 'Timestamp_ms', 'Points','Ampitude','Label'});

% - In this script, since converting depth to 3d amode data does not require 
%   estimation (just a straightforward computation), we don't need to create 
%   date/time prefix, as the result will always be the same.
% - What i meant, for depthdata, there are some process that requires
%   fine-tuning the parameters. In that case, we should create folders and
%   each folder will correspond to each of the fine-tuned parameters
% - Here, 3ddata does not require that so we can just save it, AS LONG AS 
%   we are aware to put the amode3d data to be the same directory as the 
%   corresponding depthdata, because this calculation is based on this
%   particular fine-tuned depth data
if(is_saveMat)
    % save
    mat_filename = sprintf('all_amode3d_s%s_m%s.mat', sess_str(end), meas_str);
    mat_fullpath = fullfile(path_outputdepth, mat_filename);
    save(mat_fullpath, 'all_amode3d_table');
end