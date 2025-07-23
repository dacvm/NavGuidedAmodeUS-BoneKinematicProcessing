%% HELOOOOOO

clc; clear; close all;

% [EDIT] directory to the project
path_root    = 'D:\Documents\BELANDA\PhD Thesis\Code\MATLAB\amode_navigation_experiment\experiment_b';

% [EDIT] directory to the trial
dir_trial    = "trial_0011_Session3_02";

% [EDIT] Specify folder index
folder_idx = 1;

% [EDIT] Select bone and pin
idx_bone = 2;
idx_pin  = 1;

% [EDIT]
is_display = false;

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

% Generate path to function directory
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


%% INITIALIZE SOME RIGID BODY METADATA


% 1.a) Load and organize the struct for the CT data (allBone_CT)
run('extra_structCTdata.m');

% 1.b) Get the relevant point cloud of the bone from CT scan
boneCTstl_original     = allBone_CT(idx_bone).stl;
boneCTpc_original      = pointCloud(boneCTstl_original.Points);
boneCTnormals_original = STLVertexNormals(boneCTstl_original.ConnectivityList, boneCTstl_original.Points)';
boneCTpoints_original  = [boneCTpc_original.Location'; ones(1, length(boneCTpc_original.Location))];

% 1.c) Get the relevant transformation of the bone from CT scan
T_bone_ct = allBone_CT(idx_bone).T;
bone_name = allBone_CT(idx_bone).name;
T_pin_ct  = allBone_CT(idx_bone).pin(idx_pin).T;
pin_name  = allBone_CT(idx_bone).pin(idx_pin).name;

% 2) There is this thing called REF. 
% -- It is a calibration tool for B-mode ultrasound. We were using this for 
%    bone surface reconstruction as well.
% -- They were all represented in respect of this REF. So naturally we are 
%    using this REF as our reference for everything. 
% -- Well, it is actually not really necessary to do that, we can just use 
%    Qulisys coordinate frame, but because in the data processing in the 
%    earlier experiment we are using REF coordinate frame, so let's just 
%    use it again.
ref_name   = 'B_N_REF';


%% INITIALIZE SOME DATA

% Load the rigid body data that is stored in a .csv File
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
load('all_rigidbodies_table_s3_m02_d01.mat');



%% INITIALIZE FIGURE OBJECTS AND EVERYTHING RELATED TO PLOTS

fig1 = figure('Name', 'Tibia GT in 3D');
ax1  = axes(fig1);
axis(ax1, 'equal');
hold(ax1, 'on');
grid(ax1, 'on');
xlabel(ax1, 'X');
ylabel(ax1, 'Y');
zlabel(ax1, 'Z');
view(ax1, 30,15);

fig2 = figure('Name', 'Plot Rotation Euler');
ax2  = axes(fig2);
hold(ax2, 'on');
grid(ax2, 'on');
xlabel(ax2, 'Timestamp');
ylabel(ax2, 'degree');

%% MAIN PROGRAM


% Initialize the number of timestmap
n_timestamp     = size(all_rigidbodies_table, 1);

% Allocate memory to store rotation in euler angles
rs_eul = zeros(n_timestamp, 3);

% Initialize waitbar, so the user not get bored
h = waitbar(0, 'Please wait...');

% - So, i want to get the T_boneGT_ref, that is the actual T for the bone
%   when we are doing the experiment.
% - I want to extract the most significant motion of the bone from that 
%   pose (can be flexion-extention, internal-external). 
% - I need to loop through all the T_pin_Q that has been recorded, convert
%   until i got the T_boneGT_ref, decompose to r in euler angles, and
%   collect it
for idx_t = 1:n_timestamp

    % Show the waitbar progress
    waitbar(idx_t/n_timestamp, h, sprintf('Timestamp processed: %d/%d', idx_t, n_timestamp));

    % 0) Delete every object in the plot before we go
    delete(findobj(ax1, 'Tag', 'cs_bonepin'));
    delete(findobj(ax1, 'Tag', 'cs_bonegt'));
    delete(findobj(ax1, 'Tag', 'bonesurface_ct'));

    % 1) Get the current T_pin_Q
    T_pin_Q = all_rigidbodies_table.(pin_name)(idx_t);
    T_pin_Q = T_pin_Q{1}.T;
    T_ref_Q = all_rigidbodies_table.(ref_name)(idx_t);
    T_ref_Q = T_ref_Q{1}.T;

    % 2) Calculate transformation of pin in ref 
    T_pin_ref = inv(T_ref_Q) * T_pin_Q;
    
    % 3) Calculate transformation of boneGT in ref
    T_bone_pin    = inv(T_pin_ct) * T_bone_ct;
    T_boneGT_ref  = T_pin_ref * T_bone_pin;

    % 4) Transform the point cloud
    T_ct_bone        = inv(T_bone_ct);
    boneCTpoints_GT  = T_boneGT_ref           * inv(T_bone_ct)      * boneCTpoints_original;
    boneCTnormals_GT = T_boneGT_ref(1:3, 1:3) * T_ct_bone(1:3, 1:3) * boneCTnormals_original;
    boneCTtri_GT     = triangulation(boneCTstl_original.ConnectivityList, boneCTpoints_GT(1:3,:)');

    % 5) Get the euler rotation 
    rs_eul(idx_t, :) = rad2deg(rotm2eul(T_boneGT_ref(1:3, 1:3), 'XYZ'));

    % 5) Display (This display thing is only for sanity check)
    if(is_display)
        display_axis(ax1, T_pin_ref(1:3, 4), T_pin_ref(1:3, 1:3), 30, 'T_pin_ref', 'Tag', 'cs_bonepin');
        display_axis(ax1, T_boneGT_ref(1:3, 4), T_boneGT_ref(1:3, 1:3), 30, 'T_boneGT_ref', 'Tag', 'cs_bonegt');
        trisurf(boneCTtri_GT, 'FaceColor', '#bdc3c7', 'FaceAlpha', 0.15, 'EdgeColor', 'k', 'EdgeAlpha', 0.1, 'Parent', ax1, 'Tag', 'bonesurface_ct');

        % plot the collection of r euler
        plot(ax2, rs_eul(1:idx_t,3), '-b', 'LineWidth', 2);
        drawnow;
    end

end

% Close the waitbar
close(h);

% In this particular case, we can choose which degree of freedom that we
% want to use to use as a reference to find the cycle pattern
x = 1:n_timestamp;
y = rs_eul(:,3);
plot(ax2, x, y, '-b', 'LineWidth', 2);

% Now, from this particular, y, i want to find the find the "valley" of 
% this cyclic motion. I can derivative, but why not use findpeaks function?
[pksNeg, locs] = findpeaks(-y, x, ...                    % -y turns valleys into peaks
                           'MinPeakProminence', 5, ...   % ignore tiny wiggles (tune)
                           'MinPeakDistance', 80);       % enforce separation (tune)

% Get the necessary value
valley_x = locs;          % x–coordinates of each valley
valley_y = -pksNeg;       % corresponding y–values
plot(ax2, valley_x, valley_y, 'or', 'MarkerFaceColor', 'r');

% % Write the value to a csv file   
% tmp_str = split(dir_trial, '_');
% sess_str = tmp_str{3};
% meas_str = tmp_str{4};
% csv_filename = sprintf('cycle_timestamp_s%s_m%s.csv', sess_str(end), meas_str);
% csv_fullpath = fullfile(path_outputs, csv_filename);
% writematrix(valley_x, csv_fullpath);




