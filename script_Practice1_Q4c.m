%% KYA314 - Practice 1, Q4c
% simulate the recursive sequence for varying values of lambda
clear;
close all;
clc;

% define map
f=@(x,lambda) (1+lambda).*x - x.^2;

%% Period 2 solution
% set parameters
lambda = 2.1;
Nsteps = 1000;
x0 = 0.5;

% create empty solution vector 
xtraj = NaN(Nsteps+1,1);
xtraj(1,:) = x0;

% iterate map
for i = 1:Nsteps
    x0 = f(x0,lambda);
    xtraj(i+1,:) = x0;
end

% plot solution
figure(1); clf;
plot(xtraj(100:end-1,:),xtraj(101:end,:),'.','MarkerSize',8,'Linewidth',3)
axis([-abs(lambda)-1 abs(lambda)+1 -abs(lambda)-1 abs(lambda)+1])

figure(2); clf;
plot(linspace(101,Nsteps+1,Nsteps-99),xtraj(101:end,:),'.','MarkerSize',8,'Linewidth',3)
xlim([100 1000])

%% Period 4 solution
% set parameters
lambda = 2.52;
Nsteps = 1000;
x0 = 0.5;

% create empty solution vector 
xtraj = NaN(Nsteps+1,1);
xtraj(1,:) = x0;

% iterate map
for i = 1:Nsteps
    x0 = f(x0,lambda);
    xtraj(i+1,:) = x0;
end

% plot solution
figure(1); clf;
plot(xtraj(100:end-1,:),xtraj(101:end,:),'.','MarkerSize',8,'Linewidth',3)
axis([-abs(lambda)-1 abs(lambda)+1 -abs(lambda)-1 abs(lambda)+1])

figure(2); clf;
plot(linspace(101,Nsteps+1,Nsteps-99),xtraj(101:end,:),'.','MarkerSize',8,'Linewidth',3)
xlim([100 1000])

%% Period 8 solution
% set parameters
lambda = 2.56;
Nsteps = 1000;
x0 = 0.5;

% create empty solution vector 
xtraj = NaN(Nsteps+1,1);
xtraj(1,:) = x0;

% iterate map
for i = 1:Nsteps
    x0 = f(x0,lambda);
    xtraj(i+1,:) = x0;
end

% plot solution
figure(1); clf;
plot(xtraj(100:end-1,:),xtraj(101:end,:),'.','MarkerSize',8,'Linewidth',3)
axis([-abs(lambda)-1 abs(lambda)+1 -abs(lambda)-1 abs(lambda)+1])

figure(2); clf;
plot(linspace(101,Nsteps+1,Nsteps-99),xtraj(101:end,:),'.','MarkerSize',8,'Linewidth',3)
xlim([100 1000])

%% Chaos
% set parameters
lambda = 2.7;
Nsteps = 1000;
x0 = 0.5;

% create empty solution vector 
xtraj = NaN(Nsteps+1,1);
xtraj(1,:) = x0;

% iterate map
for i = 1:Nsteps
    x0 = f(x0,lambda);
    xtraj(i+1,:) = x0;
end

% plot solution
figure(1); clf;
plot(xtraj(100:end-1,:),xtraj(101:end,:),'.','MarkerSize',8,'Linewidth',3)
axis([-abs(lambda)-1 abs(lambda)+1 -abs(lambda)-1 abs(lambda)+1])

figure(2); clf;
plot(linspace(101,Nsteps+1,Nsteps-99),xtraj(101:end,:),'.','MarkerSize',8,'Linewidth',3)
xlim([100 1000])