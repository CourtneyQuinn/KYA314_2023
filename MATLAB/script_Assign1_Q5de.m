%% KYA314 - Assignment 1, Q5de
% simulate the Henon Map for given values of alpha and beta
clear;
close all;
clc;

% define map
f=@(x,p) [1-p(1).*x(1,:).^2+x(2,:);
          p(2).*x(1,:)];

%% Period 8 solution
% set parameters
alpha = 1.25;
beta = 0.3;
p = [alpha;beta];
Nsteps = 1000;
x0 = 0.5;
y0 = 0.2;
xtemp = [x0;y0];

% create empty solution vector 
xtraj = NaN(2,Nsteps+1);
xtraj(:,1) = xtemp;

% iterate map
for i = 1:Nsteps
    xtemp = f(xtemp,p);
    xtraj(:,i+1) = xtemp;
end

% plot solution
figure(1); clf;
plot(xtraj(1,200:end),xtraj(2,200:end),'.','MarkerSize',12,'Linewidth',3)
axis([-2 2 -0.5 0.5])
set(gca,'FontSize',12)
xlabel('x_n')
ylabel('y_n')
title("Trajectories in Phase Plane")

figure(2); clf;
plot(linspace(0,Nsteps,Nsteps+1),xtraj(1,:),'.','MarkerSize',8,'Linewidth',3)
xlim([0 1000])
set(gca,'FontSize',12)
xlabel('n')
ylabel('x_n')
title("Solution Trajectories")
%% Chaos
% set parameters
alpha = 1.4;
beta = 0.3;
p = [alpha;beta];
Nsteps = 1000;
x0 = 0.5;
y0 = 0.2;
xtemp = [x0;y0];

% create empty solution vector 
xtraj = NaN(2,Nsteps+1);
xtraj(:,1) = xtemp;

% iterate map
for i = 1:Nsteps
    xtemp = f(xtemp,p);
    xtraj(:,i+1) = xtemp;
end

% plot solution
figure(1); clf;
plot(xtraj(1,200:end),xtraj(2,200:end),'.','MarkerSize',12,'Linewidth',3)
axis([-2 2 -0.5 0.5])
set(gca,'FontSize',12)
xlabel('x_n')
ylabel('y_n')
title("Trajectories in Phase Plane")

figure(2); clf;
plot(linspace(0,Nsteps,Nsteps+1),xtraj(1,:),'.','MarkerSize',8,'Linewidth',3)
xlim([0 1000])
set(gca,'FontSize',12)
xlabel('n')
ylabel('x_n')
title("Solution Trajectories")
