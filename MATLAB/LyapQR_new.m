function [ lambda,Rdiag,x,Lambda ] = LyapQR_new( M,xini,N,dM,h)
% Compute Lyapunov exponents for either a map or an array of tangent linear
% propagators
% Input
% ----------
% M : function handle or array (n,n,N+1)
%     Takes input x
% xini : array (n,)
%        initial point
% N : integer
%     number of steps over which to compute exponents
% Returns
% -------
% lams : array (n,)
%        final Lyapunov exponents
% Rdiag : array (n,N)
%         diagonal pf R at each step
% Lambda : array (n,N)
%          Lyapunov exponents at increasing windows
% x : array (n,N+1)
%     state space trajectory

n = size(xini,1);
Q = eye(n);
Rdiag = NaN(n,N);
Lambda = NaN(n,N);
x = NaN(n,N+1);

if ~exist('dM')
    hjac = 1e-5;
    dM =@(x) MyJacobian(M,x,hjac);
end


if isa(M,'function_handle')
    x(:,1) = M(xini);
    for j = 1:N
        x(:,j+1) = M(x(:,j));
        Aj = dM(x(:,j+1));
        B = Aj * Q;
        [Q,Rj] = qr(B);
        Rj = diag(Rj);
        l = find(Rj<0);
        if ~isempty(l)
            Q(:,l) = -1*Q(:,l);
        end
        Rj = abs(Rj);
        Rdiag(:,j) = Rj;
        Lambda(:,j) = 1./(j*h).*sum(log(Rdiag(:,1:j)),2);
    end
else
    for j = 1:N
        [Q,Rj] = qr(M(:,:,j+1)*Q);
        Rj = diag(Rj);
        l = find(Rj<0);
        if ~isempty(l)
            Q(:,l) = -1*Q(:,l);
        end
        Rj = abs(Rj);
        Rdiag(:,j) = Rj;
        Lambda(:,j) = 1./(j*h).*sum(log(Rdiag(:,1:j)),2);
    end
end 

lambda = Lambda(:,end);

end