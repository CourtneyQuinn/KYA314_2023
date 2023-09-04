%% Script to test MyIVP
clear;
clc;
close all;

%% x''=-x-x'+cos(t) -> solution is [sin(t);cos(t)]
% can turn above 2nd order ODE into system of two 1st order ODEs
f=@(t,x)[x(2);... 
-x(1)-x(2)+cos(t)];

%% Define MyIVP inputs
tspan=[0,2*pi]; 
h=0.05; 
xini=[0;1]; 

%% Check function with initial condition
f(0,xini)

%% Solve ODE system
[xt,t,xend]=MyIVP(f,xini,tspan,h);

%% Plot solution
plot(t,xt,'.-');
xlim([0,2*pi])
legend(['x_1';'x_2'])