function formJ(obj)
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
    [-(1+a), -(1+a), -(1-a), -(1-a);
     -(1+a), -(1+a), -(1-a), -(1-a);
     +(1-a), +(1-a), +(1+a), +(1+a);
     +(1-a), +(1-a), +(1+a), +(1+a)];
Jinv1 = obj.id;
Jinv2 = obj.id;
Jinv3 = obj.id;
Jinv4 = obj.id;
J1 = inv(Jinv1);
J2 = inv(Jinv2);
J3 = inv(Jinv3);
J4 = inv(Jinv4);
obj.J = [J1, J2, J3, J4];
end