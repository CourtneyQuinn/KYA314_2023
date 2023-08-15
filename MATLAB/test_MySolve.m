%% Script to test MySolve
clear;
clc;
format long;

%% 1D example
% Define function
f=@(x)exp(x)+x;

% Define Jacobian (numerical)
h=1e-6;
df=@(x)MyJacobian(f,x,h);

% Define MySolve inputs
x0 = 1;
tol = 1e-8;
maxit = 10;

fprintf("Numerical Jacobian \n\n")
[x,conv] = MySolve(f,x0,df,tol,maxit)

% Define Jacobian (analytical)
df=@(x)exp(x)+1;

fprintf("Analytical Jacobian \n\n")
[x,conv] = MySolve(f,x0,df,tol,maxit)


%% Rosenbrock's banana
clear;
clc;
% we find the minimum of 
f=@(x)((1-x(1))^2+100*(x(2)-x(1)^2)^2);
% by solving grad f(x)=0 with a Newton iteration:

% Define grad f(x) numerically
h=1e-5;
f2=@(x)MyJacobian(f,x,h)';   % f2 is grad f
df=@(x)MyJacobian(f2,x,h);

% Define MySolve inputs
tol=1e-6;
maxit=10;
x0=[0;0];                     % x0 is the initial guess

fprintf("Numerical Jacobians \n\n")
[x,conv]=MySolve(f2,x0,df,tol,maxit)
% the solution should be [1;1] after 5 (or so) iterations

% Define analytical Jacobian
f2 =@(x) [-2.*(1-x(1,:))-400.*x(1,:).*(x(2,:)-x(1,:).^2);... 
          200.*(x(2,:)-x(1,:).^2)];
df =@(x) [2-400.*(x(2,:)-x(1,:).^2)+800.*x(1,:).^2, -400.*x(1,:);...
          -400.*x(1,:), 200];


fprintf("Analytical Jacobians \n\n")
[x,conv]=MySolve(f2,x0,df,tol,maxit)
% the solution should be [1;1] after 5 (or so) iterations
