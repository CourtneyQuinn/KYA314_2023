function [x,converged,J] = MySolve(f,x0,df,tol,maxit)
%Solving nonlinear system of equations of form 0=f(x)
%   f has n-dimensional input and output
J = df(x0);
x1 = x0 - J\f(x0);
ek1 = norm(x1 - x0);
rk = norm(f(x0));
fprintf('Error is %g, Residual is %g\n', ek1,rk)
x0 = x1;
for i = 1:maxit
    if i == maxit
        converged = false;
        x = x0;
        break
    elseif ek1 < tol && rk < tol
        converged = true;
        J = df(x0);
        x = x0;
        break
    else
        x1 = x0 - (df(x0))\f(x0);
        ek1 = norm(x1 - x0);
        rk = norm(f(x0));
        x0 = x1;
        fprintf('Error is %g, Residual is %g\n', ek1,rk)
    end
end
if ek1 < tol && rk < tol
    converged = true;
    J = df(x0);
    x = x0;
end
end