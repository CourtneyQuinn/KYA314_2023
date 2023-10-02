%% KYA314 - Chaotic Map
% simulate the recursive sequence for varying values of r
% calculate Lyapunov exponents
clear;
close all;
clc;

%% define map
f =@(x,p) exp(-p(1).*x.^2)+p(2);

%% alpha = 5, beta in [-1,1]
% set constant parameters
alpha = 5;
Nsteps = 5000;

%% Plot Lyap exp for varying beta
% define tangent function
tan_func =@(x) -2.*x.*alpha.*exp(-alpha.*x.^2);
% take log
logdf =@(x) log(abs(tan_func(x)));

% choose beta values
beta_vals = linspace(-1,1,201);

% create empty vector for exponents
Lyaps = NaN(length(beta_vals),1);

for j = 1:length(beta_vals)
    beta = beta_vals(j);
    
    % create empty solution vector 
    xtraj = NaN(Nsteps+1,1);
    
    % input initial condition
    x0 = 0.1;
    xtraj(1,:) = x0;
    
    % iterate map
    for i = 1:Nsteps
        x0 = f(x0,[alpha;beta]);
        xtraj(i+1,:) = x0;
    end

    % We don't need to redefine as beta doesn't appear in tangent function
    Lyap = mean(logdf(xtraj));

    Lyaps(j) = Lyap;
end

figure(3); clf;
hold on;
plot([-1,1],[0,0],'k-','Linewidth',2)
plot(beta_vals,Lyaps,'r.','MarkerSize',8,'Linewidth',3)
xlabel("\beta")
ylabel("\lambda")
title("Lyapunov exponents")

%% Period 3 solution
% could not find in given alpha range - choose alpha = 7
% beta value with LE near 0 in chaotic regime
alpha = 7;
beta = -0.57;

xtraj = NaN(Nsteps+1,1);
    
% input initial condition
x0 = 0.1;
xtraj(1,:) = x0;
    
% iterate map
for i = 1:Nsteps
    x0 = f(x0,[alpha;beta]);
    xtraj(i+1,:) = x0;
end


% plot solution
figure(1); clf;
plot(xtraj(100:end-1,:),xtraj(101:end,:),'.','MarkerSize',10,'Linewidth',3)
axis([-1 1 -1 1])
xlabel("x_n")
ylabel("x_{n+1}")
title("Solution trajectory in phase space")

figure(2); clf;
plot(linspace(101,Nsteps+1,Nsteps-99),xtraj(101:end,:),'.','MarkerSize',8,'Linewidth',3)
xlim([100 1000])
xlabel("n")
ylabel("x_n")
title("Solution trajectory")