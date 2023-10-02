function dx = Stommel(t,x,p)
%Stommel Box model as defined in Dijkstra 2013
%   RHS of ODE system
dx1 = -p(2,:).*(x(1,:)-1)-x(1,:).*(1+p(3,:).^2.*(x(1,:)-x(2,:)).^2); % temperature anomalies
dx2 = p(1,:)-x(2,:).*(1+p(3,:).^2.*(x(1,:)-x(2,:)).^2); % salinity anomalies

dx = [dx1;dx2];
end