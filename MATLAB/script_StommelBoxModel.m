%% Stommel Box Model
clc;
clear;
close all;

%% Define the model
% See additional function files
%% Find equilibria
F = 0; % Fresh water flux (forcing)
alpha = 100; % timescale ratio - diffusion / temperature relaxation
mu = sqrt(6.2); % coupling parameter

p = [F;alpha;mu];

%% bifurcation diagram of unforced model
% equilibrium 1
x1 = 1;
x2 = 1;

% Set up for initial value problem solver

x0 = [x1;x2];
tspan = [0,100];
h = 0.01;

% Solve the ODE

[X,t,xeq1] = MyIVP(@(t,x)Stommel(t,x,p),x0,tspan,h);

figure(1); plot(X(1,:),X(2,:),'Linewidth',2);
set(gca,'FontSize',16)
title('Phase plane: F=0');
xlabel('x');
ylabel('y');


%% 1 Parameter continuation
RHS =@(y)Stommel(0,[y(1,:);y(2,:)],[y(3,:);alpha;mu]);

h=1e-6;
df=@(y)MyJacobian(RHS,y,h);

% Define MyTrackCurve initial points
y0 = [xeq1;F];
ytan0 = [0;0;1]; % guess in the direction of increasing F

% Solve for the true initial tangent
df0 = df(y0);
zeros0 = zeros([2 1]);
zeros1 = [zeros0;1];
f0 = [df0;transpose(ytan0)];
z = f0\zeros1;
mult = sign(transpose(z)*ytan0);
ytan = z/norm(z,Inf)*mult;

ylist1 = MyTrackCurve(RHS,df,y0,ytan,'nmax',500);

% guess in direction for decreasing F
ytan0 = [0;0;-1];

% Solve for the true initial tangent
f0 = [df0;transpose(ytan0)];
z = f0\zeros1;
mult = sign(transpose(z)*ytan0);
ytan = z/norm(z,Inf)*mult;

ylist2 = MyTrackCurve(RHS,df,y0,ytan,'nmax',500);


%% Plot tracked curve
ylist = [flip(ylist2,2),ylist1];

figure(1); clf;
hold on;
plot(ylist(3,:),ylist(2,:),'linewidth',2)
xlabel('F')
ylabel('y - salinity gradient')
set(gca,'FontSize',16)
box on;

%% Plot as bifurcation diagram
% empty vectors for eigenvalues and stability
stab = NaN(size(ylist,2),1);

for i = 1:size(ylist,2)    
    pi =[ylist(3,i);alpha;mu];
    J = MyJacobian(@(x)Stommel(0,x,pi),ylist(1:2,i),h);
    max_eig = max(eig(J(:,:,1)));
    stab(i) = max_eig > 0;
    if abs(max_eig) < 1e-2
        stab(i) = 10;
    end
end

% define dark green color for Saddle-Node bifurcation marker
dg = [77 149 66]./225;

figure(2); clf;
hold on;
for i = 1:size(ylist,2)
    if stab(i) == 0
        p1 = plot(ylist(3,i),ylist(2,i),'b.','Markersize',8);
    elseif stab(i) == 1
        p2 = plot(ylist(3,i),ylist(2,i),'r.','Markersize',8);
    elseif stab(i) == 10
        p3 = plot(ylist(3,i),ylist(2,i),'o','color',dg,'Markersize',8,'Linewidth',3);
    else
        print('Error - undefined stability type')
    end
end
xlabel('F')
ylabel('y - salinity gradient')
set(gca,'FontSize',16)
box on;
legend([p1,p2,p3],'stable eq','unstable eq','SN','location','northwest')

%% 2 Parameter Continuation
%function J = df_2param(y)
%    J = MyJacobian(@(y)Stommel(0,[y(1,:);y(2,:)],[y(3,:);y(4,:);mu]),y,h);
%end

RHS_fold =@(y) [Stommel(0,[y(1,:);y(2,:)],[y(3,:);y(4,:);mu]);
                det(df_2param(y))];

df_fold=@(y)MyJacobian(RHS_fold,y,h);

% Define MyTrackCurve initial points
f_inds = find(stab == 10);
y0 = [ylist(:,f_inds(1));alpha];
ytan0 = [0;0;0;1]; % guess i the direction of increasing params

% Solve for the true initial tangent
df0 = df_fold(y0);
zeros0 = zeros([3 1]);
zeros1 = [zeros0;1];
f0 = [df0;transpose(ytan0)];
z = f0\zeros1;
mult = sign(transpose(z)*ytan0);
ytan = z/norm(z,Inf)*mult;

yfold1 = MyTrackCurve(RHS_fold,df_fold,y0,ytan,'stepsize',1,'nmax',200);

% guess in direction for decreasing params
ytan0 = [0;0;0;-1];

% Solve for the true initial tangent
f0 = [df0;transpose(ytan0)];
z = f0\zeros1;
mult = sign(transpose(z)*ytan0);
ytan = z/norm(z,Inf)*mult;

yfold2 = MyTrackCurve(RHS_fold,df_fold,y0,ytan,'stepsize',1,'nmax',200);

%% 2 Param bifurcation diagram
yfold = [flip(yfold2,2),yfold1];

figure(3); clf;
hold on;
plot(yfold(3,:),yfold(4,:),'color',dg,'linewidth',2)
xlabel('F')
ylabel('\alpha')
set(gca,'FontSize',16)
box on;
axis([0.7 1.3 0 100])

%% add cusp point
plot(0.797,6.481,'k*','MarkerSize',8,'Linewidth',2)