clc; clear; close all

runge = @(x) 1./(1+25*x.^2);
theta = linspace(pi,0,1e4);
x = cos(theta);
y = runge(x);
k = 500;
coefficientsAndResults = polyfitOrtho(x,y,k);

intitialTime = cputime;
kk = 10:10:500-1;
mm = 10000:10000:100000;
times = nan(numel(kk),numel(mm));
for j = 1:numel(mm)
    m = mm(j);
    x = linspace(-1,1,m);
    for i = 1:numel(kk)
        k = kk(i);
        s = structfun(@(x) x(1:k),coefficientsAndResults,'UniformOutput',false);
        f = @() polyvalOrtho(x,s);
        times(i,j) = timeit(f);
    end
end
elapsedTime = cputime - intitialTime

surf(kk,mm,times')

[kkMatrix, mmMatrix] = ndgrid(kk,mm);
c = [ones(numel(kkMatrix),1) kkMatrix(:) mmMatrix(:) kkMatrix(:).*mmMatrix(:)]\times(:);
hold all
fsurf(@(k,m) c(1)+c(2)*k+c(3)*m+c(4)*k.*m,[xlim ylim])

% time complexity is O(k*m), linear in k, linear in m, bilinear when both change.