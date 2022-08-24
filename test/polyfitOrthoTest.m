clc; clear; close all


%% Runge
runge = @(x) 1./(1+25*x.^2);
m = 500;
theta = linspace(pi,0,m);
x = cos(theta); % Chebyshev points
y = runge(x);
k = m-1;
resultsTable = polyfitOrtho(x,y,k);
fplot(@(x) runge(x) - polyvalOrtho(x,resultsTable),[min(x) max(x)])

%% abs
m = 532;
theta = linspace(pi,0,m);
x = cos(theta);
y = abs(x);
k = 531;
resultsTable = polyfitOrtho(x,y,k);
figure
fplot(@(x) abs(x) - polyvalOrtho(x,resultsTable),[min(x) max(x)])

