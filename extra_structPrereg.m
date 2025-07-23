%% 2) LOAD THE PREREGISTRATION AREA

% if(is_usenavigationdata)
%     % [edit] specify your preregistration area for data with navigation
%     % (processed by GUI-based preregistration area selection script)
%     file_femurPrereg = "areasphere_femur_02-18-2025_22-14_withnav.mat";
%     file_tibiaPrereg = "areasphere_tibia_02-18-2025_22-18_withnav.mat";
% else
%     % [edit] specify your preregistration area for data without navigation
%     % (processed by GUI-based preregistration area selection script)
%     file_femurPrereg = "areasphere_femur_03-17-2025_10-40_withoutnav.mat";
%     file_tibiaPrereg = "areasphere_tibia_03-17-2025_11-49_withoutnav.mat";
% end

file_femurPrereg = "areasphere_femur_02-18-2025_22-14.mat";
% file_tibiaPrereg = "areasphere_tibia_02-18-2025_22-18.mat";
file_tibiaPrereg = "areasphere_tibia_07-23-2025_17-59.mat";

% load femur
load(fullfile(fullfile(path_bonestl, "prereg", file_femurPrereg)));

% let's structure the data, first for femur
allBone_preReg(1).name            = 'femur';
allBone_preReg(1).areas(1).name   = 'A_F_GTR';
allBone_preReg(1).areas(1).points = preregistrationArea{1}';
allBone_preReg(1).areas(1).sphere = preregistrationSphere(1,:); % xyzr
allBone_preReg(1).areas(2).name   = 'A_F_LEP';
allBone_preReg(1).areas(2).points = preregistrationArea{2}';
allBone_preReg(1).areas(2).sphere = preregistrationSphere(2,:); % xyzr
allBone_preReg(1).areas(3).name   = 'A_F_MEP';
allBone_preReg(1).areas(3).points = preregistrationArea{3}';
allBone_preReg(1).areas(3).sphere = preregistrationSphere(3,:); % xyzr

% load tibia
load(fullfile(fullfile(path_bonestl, "prereg", file_tibiaPrereg)));

% let's structure the data, second one for tibia
allBone_preReg(2).name            = 'tibia';
allBone_preReg(2).areas(1).name   = 'A_T_LEP';
allBone_preReg(2).areas(1).points = preregistrationArea{1}';
allBone_preReg(1).areas(1).sphere = preregistrationSphere(1,:); % xyzr
allBone_preReg(2).areas(2).name   = 'A_T_MEP';
allBone_preReg(2).areas(2).points = preregistrationArea{2}';
allBone_preReg(1).areas(2).sphere = preregistrationSphere(2,:); % xyzr
allBone_preReg(2).areas(3).name   = 'A_T_MAL';
allBone_preReg(2).areas(3).points = preregistrationArea{3}';
allBone_preReg(1).areas(3).sphere = preregistrationSphere(3,:); % xyzr

% clear some variables
clearvars preregistrationArea preregistrationSphere;
