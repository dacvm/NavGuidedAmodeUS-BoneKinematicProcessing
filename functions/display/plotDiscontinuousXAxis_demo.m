%% DEMO: Discontinuous x-axis with synthetic data
% This script generates two separated time windows of data and plots them
% side-by-side on one axes with a compressed gap.

clear; close all; clc;

% ----- Synthetic data (two windows separated in time) --------------------
rng(1);                     % reproducibility
fs   = 50;                  % "sample rate" (points per second)
t1   = (0 : 1/fs : 30)';    % first window:   0–30 s
t2   = (76.3 : 1/fs : 95)'; % second window:  ~76–95 s (creates the gap)

% Make a smooth oscillatory signal with small noise; similar shape for both
f0 = 0.18;                  % Hz
y1 = 13*sin(2*pi*f0*t1) - 20 + 0.6*randn(size(t1));
y2 = 13*sin(2*pi*f0*t2 + 0.35) - 20 + 0.6*randn(size(t2));

% ----- Plot using the function -------------------------------------------
figure('Color','w');
ax = axes('NextPlot','add');
title('Anterior(+) - Posterior(-)');
xlabel('Time (s)'); ylabel('mm');

% Draw with: narrow visual gap, integer tick labels, custom line style
plotDiscontinuousXAxis(gca, t1, y1, t2, y2, ...
    'GapFraction', 0.03, ...   % ~3% of total span becomes the visible gap
    'NumTicks',    10, ...     % evenly spaced ticks across the whole axis
    'TickFormat',  '%.0f', ... % integer tick labels
    'LineWidth',   1.8, ...
    'LineStyle',   '-', ...
    'Color',       [0 0.45 0.74]);  % MATLAB default blue

grid on;
