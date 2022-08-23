
runge = @(x) 1./(1+25*x.^2);
theta = linspace(pi,0,150);
x = cos(theta); % Chebyshev points
y = runge(x);
k = 100;
resultsTable = polyfitOrtho(x,y,k);

%% 