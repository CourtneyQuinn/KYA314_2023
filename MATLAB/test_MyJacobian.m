%% Script to test MyJacobian
clear;
clc;
format long;

% Define test function
func =@(x) [sin(x(3,:)) + sin(x(1,:) .* x(2,:));
    cos(x(1,:) + x(2,:) .* x(3,:).^2)];

% Define MyJacobian inputs
h = 1e-6;
x0 = [0.0,0.0,0.0];
x1 = [1.0,1.0,1.0];
xs = [x0;x1].';

if size(xs) ~= [3,2]
    fprintf(2,"Error - points array incorrectly defined\n")
end

% Calculate Jacobians
df = MyJacobian(func,xs,h);

% print individually
fprintf('df0 =\n')
disp(df(:,:,1))

fprintf('df1 =\n')
disp(df(:,:,2))
