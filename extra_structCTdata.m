%% LOAD EVERYTHING RELATED TO CT SCAN (STLS, PINS, AND CORDINATE SYSTEMS)

% [edit] specify your bone stl files 
% (processed by Mimics + 3Matic)
file_femurSTL  = "Femur_adaptiveremeshed.stl";
file_tibiaSTL  = "Tibia_adaptiveremeshed.stl";

% [edit] specify your bone pin files 
% (processed by script: extra_pinProcessStructure.m)
% file_bonePin = "output-ctpin_2025-02-21-19-08-16_partialball.mat";
file_bonePin = "output-ctpin_2025-03-17-13-26-07_fullball.mat";

%%

% load the stl file
femurSTL = stlread(fullfile(path_bonestl, file_femurSTL));
tibiaSTL = stlread(fullfile(path_bonestl, file_tibiaSTL));

% load bonePins
load(fullfile(path_outputs, file_bonePin)); 

% [edit] define the rigid body of femur and tibia 
% (processed by script from ERC project, an automatic ACS calculation)
% I got it from barabara, so ask her if you have question.
% T_femur_ct = [  -0.2209, -5.6731e-04, -0.9753, 22.4451;
%                 -0.9728, -0.0721,      0.2203, 91.6076;
%                 -0.0705,  0.9974,      0.0154, 1.2786e+03
%                    0,       0,         0,      1];
% T_tibia_ct = [ -0.4271, -0.0078, -0.9042, 15.3447;
%                -0.9029,  0.0579,  0.4260, 86.7012;
%                 0.0490,  0.9983, -0.0318, 1.2493e+03;
%                 0,       0,       0,       1];
% This one is updated version, i got the function myself
T_femur_ct = [-0.2155633328,     0.002185077,   -0.976487416,     22.46250024;
              -0.9739699967,    -0.072275141,    0.214845873,     91.56123493;
              -0.0701063115,     0.997382338,    0.017708070,   1278.59161071;
               0,                   0,                    0,       1];
T_tibia_ct = [ -0.425902209,	-0.007782197,	-0.904735732,	  15.35646317;
               -0.903441529,	 0.057810151,	 0.424795705,	  86.72583559;
                0.048997066,	 0.998297262,	-0.031652231,	1249.26505580;
                0,          	 0,         	 0,         	1];

% let's structure the data, first for femur
allBone_CT(1).name = 'femur';
allBone_CT(1).stl  = femurSTL;
allBone_CT(1).T    = T_femur_ct;
allBone_CT(1).pin(1).name           = 'P_F_PRO';
allBone_CT(1).pin(1).marker_centers = bonePins( find(strcmp({bonePins.name}, 'P_F_PRO')) ).marker_centers;
allBone_CT(1).pin(1).T              = bonePins( find(strcmp({bonePins.name}, 'P_F_PRO')) ).T;
allBone_CT(1).pin(2).name           = 'P_F_DIS';
allBone_CT(1).pin(2).marker_centers = bonePins( find(strcmp({bonePins.name}, 'P_F_DIS')) ).marker_centers;
allBone_CT(1).pin(2).T              = bonePins( find(strcmp({bonePins.name}, 'P_F_DIS')) ).T;

% second for tibia
allBone_CT(2).name = 'tibia';
allBone_CT(2).stl  = tibiaSTL;
allBone_CT(2).T    = T_tibia_ct;
allBone_CT(2).pin(1).name           = 'P_T_PRO';
allBone_CT(2).pin(1).marker_centers = bonePins( find(strcmp({bonePins.name}, 'P_T_PRO')) ).marker_centers;
allBone_CT(2).pin(1).T              = bonePins( find(strcmp({bonePins.name}, 'P_T_PRO')) ).T;
allBone_CT(2).pin(2).name           = 'P_T_DIS';
allBone_CT(2).pin(2).marker_centers = bonePins( find(strcmp({bonePins.name}, 'P_T_DIS')) ).marker_centers;
allBone_CT(2).pin(2).T              = bonePins( find(strcmp({bonePins.name}, 'P_T_DIS')) ).T;

% clear some variables
clearvars femurSTL tibiaSTL ...
          T_femur_ct T_tibia_ct ...
          bonePins;