function J = formJ(obj)

a = sqrt(3) / 3;

% notation of dNdksi
%       gp1      gp2      gp3      gp4
% N1 (-a, -a) (+a, -a) (+a, +a) (-a, +a)
% N2 (-a, -a) (+a, -a) (+a, +a) (-a, +a)
% N3 (-a, -a) (+a, -a) (+a, +a) (-a, +a)
% N4 (-a, -a) (+a, -a) (+a, +a) (-a, +a)
dNdksi1 = 0.25 * ...
    [-(1+a), -(1+a), -(1-a), -(1-a);
     +(1+a), +(1+a), +(1-a), +(1-a);
     +(1-a), +(1-a), +(1+a), +(1+a);
     -(1-a), -(1-a), -(1+a), -(1+a)];
dNdksi2 = 0.25 * ...
    [-(1+a), -(1-a), -(1-a), -(1+a);
     -(1-a), -(1+a), -(1+a), -(1-a);
     +(1-a), +(1+a), +(1+a), +(1-a);
     +(1+a), +(1-a), +(1-a), +(1+a)];

% coordinate of each node [node1, node2, node3, node4]
coordinate = reshape([obj.eNodes.coordinates], 2, []);

% dx/dksi1 @ gp1, gp2, gp3, gp4
% dy/dksi1 @ gp1, gp2, gp3, gp4
dxdksi1 = coordinate * dNdksi1;

% dx/dksi2 @ gp1, gp2, gp3, gp4
% dy/dksi2 @ gp1, gp2, gp3, gp4
dxdksi2 = coordinate * dNdksi2;

% inverse of Jacobian matrix
Jinv1 = [dxdksi1(:,1), dxdksi2(:,1)]';
Jinv2 = [dxdksi1(:,2), dxdksi2(:,2)]';
Jinv3 = [dxdksi1(:,3), dxdksi2(:,3)]';
Jinv4 = [dxdksi1(:,4), dxdksi2(:,4)]';

% Jacobian matrix
J1 = inv(Jinv1);
J2 = inv(Jinv2);
J3 = inv(Jinv3);
J4 = inv(Jinv4);

% 
J = [J1; J2; J3; J4];
end