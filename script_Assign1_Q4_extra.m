%% KYA314 - Assignment 1, Q4
% simulate the recursive sequence
clear;
close all;
clc;

% define map
f=@(x) x(1,:).^2+1;

%% Fixed points
Nsteps = 50;
x0 = [0.5-1i*sqrt(3)/2, 0.5+1i*sqrt(3)/2];

% create empty solution vector 
xtraj = NaN(Nsteps+1,size(x0,2));
xtraj(1,:) = x0;

% iterate map
for i = 1:Nsteps
    x0 = f(x0);
    xtraj(i+1,:) = x0;
end

% plot solution
figure(1); clf;
plot(real(xtraj),imag(xtraj),'.','MarkerSize',12,'Linewidth',3)
axis equal
set(gca,'FontSize',12)
xlabel('Re(x_n)')
ylabel('Im(x_n)')
title("Trajectories in Phase Space")

figure(2); clf;
plot(linspace(0,Nsteps,Nsteps+1),imag(xtraj),'.','MarkerSize',8,'Linewidth',3)
xlim([0 Nsteps])
set(gca,'FontSize',12)
xlabel('n')
ylabel('Im(x_n)')
title("Solution Trajectories")
%% Period 2 solution
% set parameters
Nsteps = 35;
x0 = [-0.5-1i*sqrt(7)/2, -0.5+1i*sqrt(7)/2];

% create empty solution vector 
xtraj = NaN(Nsteps+1,size(x0,2));
xtraj(1,:) = x0;

% iterate map
for i = 1:Nsteps
    x0 = f(x0);
    xtraj(i+1,:) = x0;
end

% plot solution
figure(1); hold on;
plot(real(xtraj),imag(xtraj),'.','MarkerSize',12,'Linewidth',3)
axis equal

figure(2); hold on;
plot(linspace(0,Nsteps,Nsteps+1),imag(xtraj),'.','MarkerSize',8,'Linewidth',3)

%% Allow for changing parameter
% define map
f=@(x,lambda) x.^2+lambda;

As = -1.5:0.01:1.5;
Bs = -1:0.01:1;
%create parameter vector
lambdas = NaN(length(As)*length(Bs),1);
i = 0;
j = 0;
for a = As
    i = i+1;
    for b = Bs
        j = j+1;
        % set parameters
        lambdas((i-1)*length(Bs)+ j) = a+1i*b;
    end
    j = 0;
end

Nsteps = 100;
x0 = zeros(length(lambdas),1);

% create empty solution vector 
xtraj = NaN(Nsteps+1,size(x0,1));
xtraj(1,:) = x0; 

% iterate map
for i = 1:Nsteps
    x0 = f(x0,lambdas);
    xtraj(i+1,:) = x0;
end

% plot values of lambda for which solutions which remain bounded
figure(2); clf; hold on;
for k = 1:length(lambdas)
    if abs(xtraj(end,k)) < 1e90
    plot(real(lambdas(k)),imag(lambdas(k)),'k.','MarkerSize',10,'Linewidth',3)
    end
end
axis equal
set(gca,'FontSize',12)
xlabel('Re(\lambda)')
ylabel('Im(\lambda)')
