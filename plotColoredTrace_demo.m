A = mod(linspace(0,360,200) + 20*randn(1,200), 360);
B = 50 + 25*sin(2*pi*(1:200)/100) + 5*randn(1,200);

% 1) new figure
fig = plotColoredTrace(A, B);

% 2) existing axes
fig2 = figure; ax2 = axes(fig2);
fig2 = plotColoredTrace(A, B, ax2);