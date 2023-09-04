%% KYA314 - Assignment 3, Q3b
% numerically estimate bifurcation value in Sal'nikov model
clear;
close all;
clc;

kappa = 0.1;

% define ODE system
f=@(t,x,p) [p(1)-kappa.*x(1,:).*exp(x(2,:));
          x(1,:).*exp(x(2,:))-x(2,:)];

%% Estimate bifurcation parameter

lam =@(mu) -0.5*(kappa*exp(mu/kappa)-mu/kappa+1);

hjac = 1e-6;
df =@(mu) MyJacobian(lam,mu,hjac);

tol = 1e-8;
maxit = 50;
[mu0, converged, ~] = MySolve(lam,0.01,df,tol,maxit)

%% ODE solver
tspan = [0,100];
hsol = 0.01;

mu_min = mu0-0.01;
mu_max = mu0+0.01;

x0 = [0.5;0.5];

[xout1,tout1,~] = MyIVP(@(t,x)f(t,x,mu_min),x0,tspan,hsol);
[xout2,tout2,~] = MyIVP(@(t,x)f(t,x,mu_max),x0,tspan,hsol);

%% Plot in phase space
figure(1); clf;
hold on;
plot(xout1(1,:),xout1(2,:),'linewidth',1.5)
plot(xout2(1,:),xout2(2,:),'linewidth',1.5)
set(gca,'FontSize',16)
legend('\mu<\mu_0','\mu>\mu_0')
xlabel('x')
ylabel('T')
box on;
