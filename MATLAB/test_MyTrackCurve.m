%% Script to test MySolve
clear;
clc;

%% 1D example
% Define function
f=@(y)y(1,:).^2+y(2,:).^2-1;

%% Define Jacobian (numerical)
h=1e-6;
df=@(y)MyJacobian(f,y,h);

% Define MyTrackCurve initial points
y0 = [0;1];
ytan = [1;0];

ylist = MyTrackCurve(f,df,y0,ytan,'stepsize',0.1);

figure(1); clf;
plot(ylist(1,:),ylist(2,:),'linewidth',2)
title('Numerical Jacobian')

%% Define Jacobian (analytical)
df=@(y)[2.*y(1,:), 2.*y(2,:)];

ylist = MyTrackCurve(f,df,y0,ytan,'stepsize',0.1);

figure(2); clf;
plot(ylist(1,:),ylist(2,:),'linewidth',2)
title('Analytical Jacobian')
