function ylist = MyTrackCurve(userf,userdf,y0,ytan,varargin)
%MyTrackCurve- track a curve of n equations for n+1 variables
%   0=f(y), f: R(n+1) to R(n)
default={'stepsize',0.01,'nmax',100,'tol',1e-5,'maxit',10};
options=MySetOptions(default,varargin);
s = options.stepsize;
nmax = options.nmax;
tol = options.tol;
maxit = options.maxit;
m = size(y0,1);
n = size(userf(y0),1);
ylist = zeros(m,nmax);
ylist(:,1) = y0;
for j = 1:nmax-1
    yj = y0 + s*ytan;
    fj =@(y)transpose(ytan)*(y - yj);
    F =@(y)[userf(y);fj(y)];
    df =@(y)[userdf(y);MyJacobian(fj,y,1e-5)];
    [yk,converged] = MySolve(F,yj,df,tol,maxit);
    while converged == 0
        s = max(s/2,1e-8);
        yj = y0 + s*ytan;
        fj =@(y)transpose(ytan)*(y - yj);
        F =@(y)[userf(y);fj(y)];
        [yk,converged] = MySolve(F,yj,df,tol,maxit);
    end
    s = min(s*2,options.stepsize);
    dfk = userdf(yk);
    zeros0 = zeros([n 1]);
    zeros1 = [zeros0;1];
    fk = [dfk;transpose(ytan)];
    z = fk\zeros1;
    mult = sign(transpose(z)*ytan);
    ytan = z/norm(z,Inf)*mult;
    ylist(:,j+1) = yk;
    y0 = yk;
end