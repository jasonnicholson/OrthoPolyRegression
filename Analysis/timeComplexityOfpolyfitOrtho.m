% Warning this takes 30 minutes or so to run.
% Takeaway: The computational complexity of polyfitOrtho is O(k*m), bilinear.

clc; clear; close all

runge = @(x) 1./(1+25*x.^2);

kk = 10:500-1;
mm = 500:500:10000;
times = nan(numel(kk),numel(mm));
for j = 1:numel(mm)
    m = mm(j);
    theta = linspace(pi,0,m);
    x = cos(theta);
    y = runge(x);
    for i = 1:numel(kk)
        k = kk(i);
        f = @() polyfitOrtho(x,y,k);
        times(i,j) = timeit(f);
    end
end

surf(kk,mm,times')

[kkMatrix, mmMatrix] = ndgrid(kk,mm);
c = [ones(numel(kkMatrix),1) kkMatrix(:) mmMatrix(:) kkMatrix(:).*mmMatrix(:)]\times(:)
hold all
fsurf(@(k,m) c(1)+c(2)k+c(3)*m+c(4)*k.*m,[xlim ylim])