
runge = @(x) 1./(1+25*x.^2);

%% Least Square Fit
m = 500;
theta = linspace(pi,0,m);
x = cos(theta); % Chebyshev points
y = runge(x);
k = m-1;
resultsTable = polyfitOrtho(x,y,k);
fplot(@(x) runge(x) - polyvalOrtho(x,resultsTable),[min(x) max(x)])