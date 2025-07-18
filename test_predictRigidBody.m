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
addpath(genpath('predictTransform'));


%% INITIALIZE SOME DATA

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

%%

Ts_boneGT_ref = all_TsReg_table.Ts_boneGT_ref;
timestamps    = all_TsReg_table.Timestamp_ms;
timestamps    = timestamps - timestamps(1);

% convert cell array to 3d matrix, each slice is the T
Ts_boneGT_ref = cat(3, Ts_boneGT_ref{:});

% get the translation parts
trans = Ts_boneGT_ref(1:3, 4, :);
trans = squeeze(trans);
trans_to_t1 = trans - trans(:,1);

% compute the velocity of each two consecutive frames
ds = sqrt( sum( diff(trans_to_t1,1,2).^2, 1 ) );  % mm\
dt = diff(timestamps');                     % ms
dv = ds./dt;                                % mm/ms
da = diff(dv);

fig1 = figure('Name', 'Test', 'Position', [50 100 1300 400]);
subplot(1,3,1);
plot(timestamps(2:end), ds); grid on; xlabel('ms'); ylabel('mm'); axis('tight'); hold('on');
subplot(1,3,2);
plot(timestamps(2:end), dv); grid on; xlabel('ms'); ylabel('mm/ms'); axis('tight');
subplot(1,3,3);
plot(timestamps(3:end), da); grid on; xlabel('ms'); ylabel('mm/ms^2'); axis('tight');

figure;
plot3(trans(1,:), trans(2,:), trans(3,:), '-ob'); grid on; axis equal; hold on;


% %% MAIN PROGRAM
% 
% Ts_boneGT_ref = all_TsReg_table.Ts_boneGT_ref;
% timestamps    = all_TsReg_table.Timestamp_ms;
% 
% onlydiff_error  = [];
% predicted_error = []; 
% 
% for i=5:length(Ts_boneGT_ref)-1
% 
%     tmp = Ts_boneGT_ref(i-4:i);
%     current_Tseries = cat(3, tmp{:});
%     current_tseries = timestamps(i-4:i);
% 
%     % T_predicted = predictNextTransform_v8(current_Tseries, current_tseries', 2, 2.0);
%     T_predicted = predictNextTransform_v8(current_Tseries, current_tseries');
% 
%     % Calculate without prediction
%     error_matrix = inv(Ts_boneGT_ref{i+1}) * Ts_boneGT_ref{i};
%     translation_error = norm(error_matrix(1:3,4));
%     rotation_error_rad = acos((trace(error_matrix(1:3,1:3)) - 1) / 2);
%     onlydiff_error = [translation_error, rad2deg(rotation_error_rad)];
%     % fprintf('Diff\t=[%.4f, %4f]\n', onlydiff_error(1), onlydiff_error(2));
% 
% 
%     % Calculate the error
%     error_matrix = inv(Ts_boneGT_ref{i+1}) * T_predicted;
%     translation_error = norm(error_matrix(1:3,4));
%     rotation_error_rad = acos((trace(error_matrix(1:3,1:3)) - 1) / 2);
%     predicted_error = [translation_error, rad2deg(rotation_error_rad)];
%     % fprintf('Pred\t=[%.4f, %4f]\n\n', predicted_error(1), predicted_error(2));
% 
% 
%     plot3(T_predicted(1,4), T_predicted(2,4), T_predicted(3,4), 'or'); grid on; axis equal; hold on;
% 
% end
% 
% mean(onlydiff_error, 1);
% mean(predicted_error, 1);


% for i=2:length(Ts_boneGT_ref)-1
%     T_k = Ts_boneGT_ref{i};
%     T_km1 = Ts_boneGT_ref{i-1};
% 
%     % Prediction (constant‑velocity)
%     Xi_k   = tform2se3(T_k);         % current 6×1
%     Xi_km1 = tform2se3(T_km1);       % previous 6×1
%     v      = Xi_k - Xi_km1;          % velocity in R^6
%     Xi_pred = Xi_k + v;              % next‑frame prediction
%     T_pred  = se32tform(Xi_pred);    % back → rigid3d
% 
%     % Calculate without prediction
%     error_matrix = inv(Ts_boneGT_ref{i+1}) * Ts_boneGT_ref{i};
%     translation_error = norm(error_matrix(1:3,4));
%     rotation_error_rad = acos((trace(error_matrix(1:3,1:3)) - 1) / 2);
%     onlydiff_error = [translation_error, rad2deg(rotation_error_rad)];
%     % fprintf('Diff\t=[%.4f, %4f]\n', onlydiff_error(1), onlydiff_error(2));
% 
% 
%     % Calculate the error
%     error_matrix = inv(Ts_boneGT_ref{i+1}) * T_pred;
%     translation_error = norm(error_matrix(1:3,4));
%     rotation_error_rad = acos((trace(error_matrix(1:3,1:3)) - 1) / 2);
%     predicted_error = [translation_error, rad2deg(rotation_error_rad)];
%     % fprintf('Pred\t=[%.4f, %4f]\n\n', predicted_error(1), predicted_error(2));
% end
% 
% mean(predicted_error, 1)
% 
% 
% % Convert rigid3d to 6‑vector (se(3) twist)
% function xi = tform2se3(T)
%     % T is a rigid3d object
%     R = T(1:3, 1:3);           % 3×3
%     t = T(1:3, 4)';       % 1×3 → row
%     % axis‑angle [ax ay az theta]
%     axang = rotm2axang(R);    
%     axis = axang(1:3);
%     theta = axang(4);
%     omega = axis * theta;     % 1×3
%     xi = [omega, t];          % 1×6: [ωx ωy ωz tx ty tz]
% end
% 
% % Convert 6‑vector back to rigid3d
% function T = se32tform(xi)
%     omega = xi(1:3);
%     theta = norm(omega) + eps;
%     axis = omega / theta;
%     R = axang2rotm([axis, theta]);
%     t = xi(4:6);
%     T = [R, t'; 0 0 0 1];
% end