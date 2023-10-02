function J = df_2param(y)
    mu = sqrt(6.2);
    h = 1e-6;
    J = MyJacobian(@(x)Stommel(0,x,[y(3,:);y(4,:);mu]),[y(1,:);y(2,:)],h);
    J = J(:,:,1);
end