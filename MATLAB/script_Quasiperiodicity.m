%% Forced Van der Pol Oscillator
clc;
clear;
close all;

%% Define the model
% Here we model the forced Van der Pol - Duffing oscillator
Van_Duff=@(t,x,p)[x(2,:); % velocity
            p(1,:).*(1-x(1,:).^2).*x(2,:)-p(2,:).*x(1,:)-p(3,:).*x(1,:).^3+p(4,:).*cos(x(3,:)); % main Duffing equation
            p(5,:)*x(3,:).^0]; % new variable for forcing period

%% Set parameters
% µ = 0.1,  α = 1,     β = 0, K = 0.5,  ω = 1         P
% µ = 0.1,  α = 1,     β = 0, K = 0.5,  ω = sqrt(2)   QP
% µ = 0.1,  α = 1,     β = 1, K = 0.5,  ω = 1         QP
% µ = 1.75, α = 1,     β = 0, K = 0.5,  ω = pi/5      QP
% µ = 1,    α = −1.44, β = 1, K = 0.45, ω = 1.2       Chaos

mu = 1; % nonlinear damping
alpha = -1.44; % stiffness
beta = 1; % nonlinearity
K = 0.45; % forcing strength
omega = 1.2; % frequency

p = [mu;alpha;beta;K;omega];

%% Solve IVP

x1 = 0;
x2 = 2;
x3 = 0;

% Set up for initial value problem solver

x0 = [x1;x2;x3];
tspan = [0,5000];
timescale = 2/omega*pi;
h = 0.01*timescale;

% Solve the ODE

[X,t,xeq1] = MyIVP(@(t,x)Van_Duff(t,x,p),x0,tspan,h);

figure(1); 
plot(t,X(1,:));
set(gca,'FontSize',16)
xlabel('t');
ylabel('x');

%% Plot poincare map

figure(2); 
plot(X(1,1000:100:end),X(2,1000:100:end),'.','MarkerSize',8,'Linewidth',2);
set(gca,'FontSize',16)
xlabel('x');
ylabel('v');
%axis([-1 1 0 4])

%% Plot 3D trajectory

figure(3);
plot3(X(1,1000:end),X(2,1000:end),cos(X(3,1000:end)),'.','MarkerSize',8,'Linewidth',2);
set(gca,'FontSize',16)
xlabel('x');
ylabel('v');
zlabel('z');

%% Calculate Lyapunov exponents
hjac = 1e-6;
Js = MyJacobian(@(x)Van_Duff(0,x,p),X,hjac);
M = NaN(size(Js));
for j = 1:size(Js,3)
    M(:,:,j) = expm(Js(:,:,j)*h);
end

N = size(X,2)-1;

[ lambda,~,~,~] = LyapQR_new(M,x0,N,[],h);
disp(lambda)

%% Calculate frequency spectrum

fs = 1/h;
N = length(X(1,:));
xdft = fft(X(1,:));
xdft = xdft(1:N/2+1);
psdx = (1/(fs*N)) * abs(xdft).^2;
psdx(2:end-1) = 2*psdx(2:end-1);
freq = 0:fs/length(X(1,:)):fs/2;

figure(4);
plot(timescale*freq,psdx./max(psdx))
grid on
title("Spectrum")
xlabel("Frequency")
ylabel("Power spectral density")
xlim([0 5])
