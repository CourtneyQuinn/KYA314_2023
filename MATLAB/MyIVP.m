function [xt,t,xend] = MyIVP(f,x0,tspan,h)
% Solving initial-value problems for ODEs of form x*(t)=f(t,x(t))
%   1D time input and n-D x input to n-D x output
%   4th-order Runge-Kutta

N = round((tspan(2)-tspan(1))/h);
n = size(x0,1);
m = size(x0,2);
t = NaN(N+1,1);
t(1) = tspan(1);

if m > 1
    xt = NaN(n,m,N+1);
    xt(:,:,1) = x0;
else
    xt(:,1) = x0;
end
    
for k = 1:N
    t0 = t(k);
    if m > 1
        x0 = xt(:,:,k);
    else
        x0 = xt(:,k);
    end
    j1 = f(t0, x0);
    j2 = f(t0 + (h / 2), x0 + (h / 2).*j1);
    j3 = f(t0 + (h / 2), x0 + (h / 2).*j2);
    j4 = f(t0 + h, x0 + h.*j3);
    if m > 1
        xt(:,:,k+1) = x0 + (h / 6).*(j1 + 2*j2 + 2*j3 + j4);
    else
        xt(:,k+1) = x0 + (h / 6).*(j1 + 2*j2 + 2*j3 + j4);
    end 
    t0 = tspan(1)+ k*h;
    t(k+1) = t0;
end
if nargout > 2
    if m > 1
        xend = xt(:,:,end);
    else
        xend = xt(:,end);
    end 
   
end
end