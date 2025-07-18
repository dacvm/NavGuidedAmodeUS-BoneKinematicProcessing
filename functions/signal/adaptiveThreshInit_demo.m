clear; clc; close all;

filename = 'ust_reqdata_24.mat';
load(filename);

% get the selected peak (from user)
init_depth_mm    = ust_datapacket.init;
init_depth_idx   = round(init_depth_mm/ust_usspec.ds);
init_window_mm   = 1.5;
init_window_idx  = round(init_window_mm / ust_usspec.ds);

[peak_threshold, prom_threshold] = adaptiveThreshInit(ust_datapacket.mmode_proc, init_depth_idx, init_window_idx, true);

