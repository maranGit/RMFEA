function Calculate_ElementStiffness_Force(obj, mat)
% updated Lagrange element subroutine for finite strain
% update element stiffness and Fint for current displacement
% Input: material subroutine, X, displacement
% Output: (1) Fint; (2) element stiffness
%
% local derivative of shape function
% dN1/dksi1 dN2/dksi1 dN3/dksi1 dN4/dksi1 
% dN1/dksi2 dN2/dksi2 dN3/dksi2 dN4/dksi2 
%
numgp = 4;
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
dNdksiALL = [dNdksi1; dNdksi2; dNdksi3; dNdksi4];
%
% X1, X2, X3, X4
% Y2, Y3, Y4, Y4
X = reshape([obj.eNodes.coordinates], 2, []);
% u1, v1, u2, v2, u3, v3, u4, v4
disp = [obj.nedof.v]'; % should change to edofs in future version @@@
numDof = length(disp);
% x1, x2, x3, x4
% y2, y3, y4, y4
x = X + reshape(disp, 2, []);
%
% loop over integration point for Fint and K
Fint = zeros(numDof, 1);
K = zeros(numDof, numDof);
for gp = 1:numgp
    dNdksi = dNdksiALL([2*gp-1, 2*gp], :);
    % form Jacobian matrix
    Jcurr = dNdksi * transpose(x); % Jacobian (dx/dksi)'
    dNdx = Jcurr \ dNdksi;
    Jref = dNdksi * transpose(X); % Jacobian (dX/dksi)'
    dNdX = Jref \ dNdksi;
    F = x * transpose(dNdX); % deformation gradient
    % construct a 4*8 B matrix for finite deformation
    B = zeros(4, 8);
    B(1:2, 1:2:7) = dNdx;
    B(3:4, 2:2:8) = dNdx;
    b = F * transpose(F); % left Cauchy-Green strain tensor
    b_voigt = [b(1,1); b(2,2); b(1,2) + b(2,1)];
    [C, sigma] = mat.calcStressTangent(b_voigt);
    obj.stress(:, gp) = sigma;
    obj.strain(:, gp) = 0.5 * ( b_voigt - [1; 1; 0] );
    % form geometric stiffness term
    sigma_tensor = [sigma(1), sigma(3); sigma(3), sigma(2)];
    sigma_matrix = zeros(4, 4);
    sigma_matrix(1:2, 1:2) = sigma_tensor;
    sigma_matrix(3:4, 3:4) = sigma_tensor;
    % form material stiffness term
    temp = [1, 0, 0, 0; 0, 0, 0, 1; 0, 1, 1, 0];
    C_matrix = transpose(temp) * C * temp;
    % update Fint and K
    K = K + transpose(B) * (C_matrix + sigma_matrix) * B * det(Jcurr);
    Fint = Fint + transpose(B) * [sigma(1); sigma(3); sigma(3); sigma(2)] * det(Jcurr);
end
obj.ke = K;
obj.Fint = Fint;
%
% construct element force vector (including body force, surface loading and displacement force)
% pDisp = transpose(disp); % collect displacement of all dof
% pDisp(obj.dofMap>0) = 0; % set displacement of free dof to 0
% obj.fde = K * transpose(pDisp);
obj.fde = zeros(8, 1);
obj.fee = obj.foe - obj.fde;
%
end