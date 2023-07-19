function [xout] = LogisticMap(x,lambda)
%LogisticMap 1D map exhibiting complex behaviour
% Input
% ----------
% x : array (1, N)
%     Points to iterate
%     If N = 1, x is one point
% lambda : float
%     parameter value
% Returns
% -------
% xout : array, shape (n, N)
%     Array containing 1-step interate of map
%     2nd dimension denotes iterates for different starting points
xout = lambda.*x(1,:).*(1-x(1,:));
end