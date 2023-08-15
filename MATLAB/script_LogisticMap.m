%% KYA314 - Logistic Map
% simulate the recursive sequence for varying values of lambda
clear;
close all;
clc;

% set parameters
lambda = 1.1;
x0 = 0.1;
Nsteps = 1000;

% create empty solution vector 
xtraj = NaN(Nsteps+1,1);

% input initial condition
xtraj(1,:) = x0;

% iterate map
for i = 1:Nsteps
    x0 = LogisticMap(x0,lambda);
    xtraj(i+1,:) = x0;
end

% plot solution
figure(1); clf;
plot(xtraj(100:end-1,:),xtraj(101:end,:),'.','MarkerSize',10,'Linewidth',3)
axis([-abs(lambda)-1 abs(lambda)+1 -abs(lambda)-1 abs(lambda)+1])
xlabel("x_n")
ylabel("x_{n+1}")
title("Solution trajectory in phase space")

figure(2); clf;
plot(linspace(101,Nsteps+1,Nsteps-99),xtraj(101:end,:),'.','MarkerSize',8,'Linewidth',3)
xlim([100 1000])
xlabel("n")
ylabel("x_n")
title("Solution trajectory")
