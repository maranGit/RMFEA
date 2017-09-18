function Calculate_ElementStiffness_Force(obj, mat)
% update element stiffness and Fint for current displacement
% Input: nodal displacement
% Output: (1) residual; (2) element stiffness
%{
% delete test case with linear elastic material
obj.setInternalMaterialProperties;

E1 = obj.E1;
E2 = obj.E2;
G = obj.G;

% material stiffness
D = [E1, E2, 0;
     E2, E1, 0;
     0,  0,  G];
%}
%
% local derivative of shape function
% dN1/dksi1 dN2/dksi1 dN3/dksi1 dN4/dksi1 
% dN1/dksi2 dN2/dksi2 dN3/dksi2 dN4/dksi2 
a = sqrt(3) / 3;
dNdksi1 = 0.25 * ... % @ gp1
    [-(1+a), +(1+a), +(1-a), -(1-a);
     -(1+a), -(1-a), +(1-a), +(1+a)];
dNdksi2 = 0.25 * ... % @ gp2
    [-(1+a), +(1+a), +(1-a), -(1-a);
     -(1-a), -(1+a), +(1+a), +(1-a)];
dNdksi3 = 0.25 * ... % @ gp3
    [-(1-a), +(1-a), +(1+a), -(1+a);
     -(1-a), -(1+a), +(1+a), +(1-a)];
dNdksi4 = 0.25 * ... % @ gp4
    [-(1-a), +(1-a), +(1+a), -(1+a);
     -(1+a), -(1-a), +(1-a), +(1+a)];
%
% global derivative of shape function
% dN1/dx dN2/dx dN3/dx dN4/dx 
% dN1/dy dN2/dy dN3/dy dN4/dy 
J = obj.formJ;
J1 = J(1:2, :);
J2 = J(3:4, :);
J3 = J(5:6, :);
J4 = J(7:8, :);
dNdx1 = J1 * dNdksi1;
dNdx2 = J2 * dNdksi2;
dNdx3 = J3 * dNdksi3;
dNdx4 = J4 * dNdksi4;
%
% construct B matrix
Bmat1 = [dNdx1(1, 1), 0,           dNdx1(1, 2), 0,           dNdx1(1, 3), 0,           dNdx1(1, 4), 0;
         0,           dNdx1(2, 1), 0,           dNdx1(2, 2), 0,           dNdx1(2, 3), 0,           dNdx1(2, 4);
         dNdx1(2, 1), dNdx1(1, 1), dNdx1(2, 2), dNdx1(1, 2), dNdx1(2, 3), dNdx1(1, 3), dNdx1(2, 4), dNdx1(1, 4)];
Bmat2 = [dNdx2(1, 1), 0,           dNdx2(1, 2), 0,           dNdx2(1, 3), 0,           dNdx2(1, 4), 0;
         0,           dNdx2(2, 1), 0,           dNdx2(2, 2), 0,           dNdx2(2, 3), 0,           dNdx2(2, 4);
         dNdx2(2, 1), dNdx2(1, 1), dNdx2(2, 2), dNdx2(1, 2), dNdx2(2, 3), dNdx2(1, 3), dNdx2(2, 4), dNdx2(1, 4)];
Bmat3 = [dNdx3(1, 1), 0,           dNdx3(1, 2), 0,           dNdx3(1, 3), 0,           dNdx3(1, 4), 0;
         0,           dNdx3(2, 1), 0,           dNdx3(2, 2), 0,           dNdx3(2, 3), 0,           dNdx3(2, 4);
         dNdx3(2, 1), dNdx3(1, 1), dNdx3(2, 2), dNdx3(1, 2), dNdx3(2, 3), dNdx3(1, 3), dNdx3(2, 4), dNdx3(1, 4)];
Bmat4 = [dNdx4(1, 1), 0,           dNdx4(1, 2), 0,           dNdx4(1, 3), 0,           dNdx4(1, 4), 0;
         0,           dNdx4(2, 1), 0,           dNdx4(2, 2), 0,           dNdx4(2, 3), 0,           dNdx4(2, 4);
         dNdx4(2, 1), dNdx4(1, 1), dNdx4(2, 2), dNdx4(1, 2), dNdx4(2, 3), dNdx4(1, 3), dNdx4(2, 4), dNdx4(1, 4)];
%
% update element strain and stress for each integration point
%   eps_xx_1,   eps_xx_2,   eps_xx_3,   eps_xx_4
%   eps_yy_1,   eps_yy_2,   eps_yy_3,   eps_yy_4
% 2*eps_xx_1, 2*eps_xx_2, 2*eps_xx_3, 2*eps_xx_4
%
disp = [obj.nedof.v]'; % should change to edofs in future version @@@
obj.strain = [Bmat1*disp, Bmat2*disp, Bmat3*disp, Bmat4*disp];
[D1, sigma1] = mat.calcStressTangent(obj.strain(:,1));
[D2, sigma2] = mat.calcStressTangent(obj.strain(:,2));
[D3, sigma3] = mat.calcStressTangent(obj.strain(:,3));
[D4, sigma4] = mat.calcStressTangent(obj.strain(:,4));
obj.stress = [sigma1, sigma2, sigma3, sigma4];
%
% construct element stiffness matrix
K = transpose(Bmat1) * D1 * Bmat1 / det(J1) ...
  + transpose(Bmat2) * D2 * Bmat2 / det(J2) ...
  + transpose(Bmat3) * D3 * Bmat3 / det(J3) ...
  + transpose(Bmat4) * D4 * Bmat4 / det(J4);
obj.ke = K;
obj.Bmat_1 = Bmat1;
obj.Bmat_2 = Bmat2;
obj.Bmat_3 = Bmat3;
obj.Bmat_4 = Bmat4;
%
% construct internal force
obj.Fint = transpose(Bmat1) * sigma1 / det(J1) ...
    + transpose(Bmat2) * sigma2 / det(J2) ...
    + transpose(Bmat3) * sigma3 / det(J3) ...
    + transpose(Bmat4) * sigma4 / det(J4);
%
% construct element force vector (including body force, surface loading and displacement force)
pDisp = [obj.nedof.v]; % collecte displacement of all dof
pDisp(obj.dofMap>0) = 0; % set displacement of free dof to 0
obj.fde = K * transpose(pDisp);
obj.fee = obj.foe - obj.fde;
%
end