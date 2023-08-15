%% KYA314 - Assignment 2, Q5c
% simulate the predator-prey model
clear;
close all;
clc;

% define ODE system
f=@(t,x,p) [p(1).*x(1,:)-p(2).*x(1,:).*x(2,:);
          -p(3).*x(2,:)+p(4).*x(1,:).*x(2,:)];

%% Choose parameters and IC
r = 2;
m = 2.5;
a = 1;
b = 3;

p = [r;m;a;b];

x01 = [10;10];
x02 = [1;5];
x03 = [3;1];

%% ODE solver
tspan = [0,100];

[tout1,xout1] = ode45(@(t,x)f(t,x,p),tspan,x01);
[tout2,xout2] = ode45(@(t,x)f(t,x,p),tspan,x02);
[tout3,xout3] = ode45(@(t,x)f(t,x,p),tspan,x03);

%% Plot in phase space
figure(1); clf;
hold on;
plot(xout1(:,1),xout1(:,2))
plot(xout2(:,1),xout2(:,2))
plot(xout3(:,1),xout3(:,2))
set(gca,'FontSize',16)
xlabel('G')
ylabel('H')
box on;

%% Plot in time
figure(2); clf;
hold on;
plot(tout1,xout1(:,1),'LineWidth',2)
plot(tout1,xout1(:,2),'LineWidth',2)
set(gca,'FontSize',16)
xlabel('time')
ylabel('population')
legend('Grass','Horses')
box on;

%% Interpretion
% When starting with initial conditions of 10 grass, 10 horses, we see the
% horses initially boom in population, while the grass is consumed
% entirely. With no grass to eat, the horses slowly die off.  Once the
% horses drop below a critical level, the grass can begin to grow again and
% has a sharp increase in population. Once there is a sufficient amount of
% grass the horses can begin to grow in numbers again. Thus the cycle
% restarts.