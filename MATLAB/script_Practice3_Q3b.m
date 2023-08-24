%% KYA314 - Practice 3, Q3b
% simulate the Vander Pol oscillator
clear;
close all;
clc;

% define ODE system
f=@(t,x,p) [x(2,:);
          p(1).*(1-x(1,:).^2).*x(2,:)-x(1,:)];

%% alpha < 0
alpha = -1;

x0 = [0.1;0.1];

% ODE solver
tspan = [0,100];

[~,xout1] = ode45(@(t,x)f(t,x,alpha),tspan,x0);

% Plot in phase space
figure(1); clf;
hold on;
plot(xout1(:,1),xout1(:,2))

set(gca,'FontSize',16)
xlabel('$x$','interpreter','latex')
ylabel('$\dot{x}$','interpreter','latex')
box on;

%% alpha = 0
alpha = 0;

x0 = [0.1;0.1];

% ODE solver
tspan = [0,100];

[~,xout2] = ode45(@(t,x)f(t,x,alpha),tspan,x0);

% Plot in phase space
figure(1); clf;
hold on;
plot(xout2(:,1),xout2(:,2))

set(gca,'FontSize',16)
xlabel('$x$','interpreter','latex')
ylabel('$\dot{x}$','interpreter','latex')
box on;

%% alpha > 0
alpha = 1;

x0 = [0.1;0.1];

% ODE solver
tspan = [0,100];

[~,xout3] = ode45(@(t,x)f(t,x,alpha),tspan,x0);

% Plot in phase space
figure(1); clf;
hold on;
plot(xout3(:,1),xout3(:,2))

set(gca,'FontSize',16)
xlabel('$x$','interpreter','latex')
ylabel('$\dot{x}$','interpreter','latex')
box on;

%% alpha >> 0
alpha = 10;

x0 = [0.1;0.1];

% ODE solver
tspan = [0,100];

[~,xout4] = ode45(@(t,x)f(t,x,alpha),tspan,x0);

% Plot in phase space
figure(1); clf;
hold on;
plot(xout4(:,1),xout4(:,2))

set(gca,'FontSize',16)
xlabel('$x$','interpreter','latex')
ylabel('$\dot{x}$','interpreter','latex')
box on;
