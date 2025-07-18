clear; clc; close all;

filename = 'D:\Documents\BELANDA\PhD Thesis\Code\MATLAB\amode_navigation_experiment\experiment_b\outputs\output_allestdepths\depthdata_original\ust_reqdata_25.mat';
load(filename);

%%

% get the selected peak (from user)
init_depth_mm    = 13.25968;
init_depth_idx   = round(init_depth_mm/ust_usspec.ds);
init_window_mm   = 1.5;
init_window_idx  = round(init_window_mm / ust_usspec.ds);
[peak_threshold, prom_threshold] = adaptiveThreshInit(ust_datapacket.mmode_proc, init_depth_idx, init_window_idx, true);


%%

% get the length data
length_data     = size(ust_datapacket.mmode_proc, 2);

% prepare the figure
fig1 = figure;
ax1  = axes(fig1);
imagesc(ax1, 1:length_data, ust_usspec.s_vector, ust_datapacket.mmode_proc); hold on;
ylim(ax1, [0 20]);
axis(ax1, 'xy'); 

% [7] -----------------------------------------------------------------
% Setting up all the required parameter for depth detection
% some required parameters for mmodeDepthDetection
new_seed   = init_depth_mm;
depth_unit = ust_usspec.ds;

% peak threshold parameters
peakparam.peak_threshold = peak_threshold;
peakparam.prom_threhsold = prom_threshold;

% moving window parameters
movwinparam.length_window = 500;
movwinparam.length_stride = 100;

% bone cluster parameters
clusterparam.epsilon    = 0.16;
clusterparam.minpts     = 35;
clusterparam.cov        = [300 0; 0 1.5];
outlierparam.gradthresh = 3;
outlierparam.minpoint   = 10;
smoothingparam          = 0.0001;

% detect the bone from mmode image
estimatedbone_proc = mmodeDepthDetection( ust_datapacket.mmode_proc, depth_unit, new_seed, peakparam, ...
                                          movwinparam, clusterparam, outlierparam, smoothingparam, ...
                                          ax1, true);