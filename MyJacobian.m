function df = MyJacobian(f,x,h)
% Jacobian of a function with arbitrary dimensions
% 2nd Order Central Difference
% Parameters
% ----------
% f : function handle
%     Takes input x
% x : array (n, N)
%     Points at which Jacobian is to be calculated
%     If N = 1, x is one point
% h : float
%     Step size for finite differences
% Returns
% -------
% df : array, shape (m, n, N)
%     Array containing Jacobian matrix in first 2 dimensions
%     3rd dimension denotes which point in original x array
n = size(x,1);
N = size(x,2);

m = size(f(x),1);
df = NaN(m,n,N);
for i = 1:n
    xi1 = x;
    xi2 = x;
    xi1(i,:) = x(i,:) + h;
    xi2(i,:) = x(i,:) - h;
    dfi = (f(xi1) - f(xi2)) / (2*h);
    df(:,i,:) = dfi;
end
end