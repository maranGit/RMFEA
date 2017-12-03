function Calculate_ElementStiffness_Force(obj, mat)
% Ran Ma
% 11/22/2017
% generalized element subroutine for all finite strain plastic model
% element subroutine is based on warp3d
% 
% Input: X, u_n, u_np1
% 
% pass to material subroutine: F_n, F_np1, tau_n, hardening_n
% take from material subroutine: tau_np1, hardening_np1
% 
% Output: (1) Fint; (2) element stiffness
%
% local derivative of shape function
% dN1/dksi1 dN2/dksi1 dN3/dksi1 dN4/dksi1 
% dN1/dksi2 dN2/dksi2 dN3/dksi2 dN4/dksi2 
%
% dim = 2;
numgp = 4;
sqrt3 = sqrt(3) / 3;
dNdksi1 = 0.25 * ... % @ gp1
    [-(1+sqrt3), +(1+sqrt3), +(1-sqrt3), -(1-sqrt3);
     -(1+sqrt3), -(1-sqrt3), +(1-sqrt3), +(1+sqrt3)];
dNdksi2 = 0.25 * ... % @ gp2
    [-(1+sqrt3), +(1+sqrt3), +(1-sqrt3), -(1-sqrt3);
     -(1-sqrt3), -(1+sqrt3), +(1+sqrt3), +(1-sqrt3)];
dNdksi3 = 0.25 * ... % @ gp3
    [-(1-sqrt3), +(1-sqrt3), +(1+sqrt3), -(1+sqrt3);
     -(1-sqrt3), -(1+sqrt3), +(1+sqrt3), +(1-sqrt3)];
dNdksi4 = 0.25 * ... % @ gp4
    [-(1-sqrt3), +(1-sqrt3), +(1+sqrt3), -(1+sqrt3);
     -(1+sqrt3), -(1-sqrt3), +(1-sqrt3), +(1+sqrt3)];
dNdksiALL = [dNdksi1; dNdksi2; dNdksi3; dNdksi4];
%
% X1, X2, X3, X4
% Y2, Y3, Y4, Y4
X = reshape([obj.eNodes.coordinates], 2, []);
% u1, v1, u2, v2, u3, v3, u4, v4
disp_n = [obj.nedof.vn]'; % should change to edofs in future version @@@
disp_np1 = [obj.nedof.v]';
numDof = length(disp_n);
% x1, x2, x3, x4
% y2, y3, y4, y4
x_n = X + reshape(disp_n, 2, []);
x_np1 = X + reshape(disp_np1, 2, []);
%
% loop over integration point for Fint and K
F_n = eye(3); % initialize F_n in 3d
F_np1 = eye(3); % initialize F_np1 in 3d
Fint = zeros(numDof, 1);
K = zeros(numDof, numDof);
for gp = 1:numgp
    %
    % compute shape function gradient
    dNdksi = dNdksiALL([2*gp-1, 2*gp], :);
    Jcurr = dNdksi * transpose(x_np1); % Jacobian (dx/dksi)'
    dNdx = Jcurr \ dNdksi;
    Jref = dNdksi * transpose(X); % Jacobian (dX/dksi)'
    dNdX = Jref \ dNdksi;
    % 
    % compute F_np1
    F_n(1:2, 1:2) = x_n * transpose(dNdX); % deformation gradient
    F_np1(1:2, 1:2) = x_np1 * transpose(dNdX);
    %
    % compute deformation gradient for intermediate configuration
    alp = obj.alpha;
    F_npa = alp * F_np1 + (1 - alp) * F_n;
    f_np1 = F_np1 / F_n;
    f_npa = (1 - alp) * eye(3) + alp * f_np1;
    f_npa_tilde = f_np1 / f_npa;
    %
    % compute strain rate tensor in rotated configuration: D*delta(t)
    temp = eye(3) - inv(f_np1 * transpose(f_np1));
    eps_tilde = 0.5 * transpose(f_npa_tilde) * (temp) * f_npa_tilde;
    R_npa = polarDecomp(F_npa);
    D_tensor = transpose(R_npa) * eps_tilde * R_npa;
    voigt = [1;5;9;8;7;4];
    D = D_tensor(voigt);
    D(4:6) = D(4:6) * 2;
    %{
    % another way of computing strain
    Jhalf = dNdksi * transpose((x_np1 + x_n)/2);
    dNdx1 = Jhalf \ dNdksi;
    Bhalf = [dNdx1(1, 1), 0,           dNdx1(1, 2), 0,           dNdx1(1, 3), 0,           dNdx1(1, 4), 0;
         0,           dNdx1(2, 1), 0,           dNdx1(2, 2), 0,           dNdx1(2, 3), 0,           dNdx1(2, 4);
         dNdx1(2, 1), dNdx1(1, 1), dNdx1(2, 2), dNdx1(1, 2), dNdx1(2, 3), dNdx1(1, 3), dNdx1(2, 4), dNdx1(1, 4)];
    eps = Bhalf * (disp_np1 - disp_n);
    eps_tilde = [eps(1),eps(3)/2,0;eps(3)/2,eps(2),0;0,0,0];
    R_npa = polarDecomp(F_npa);
    D_tensor = transpose(R_npa) * eps_tilde * R_npa;
    voigt = [1;5;9;8;7;4];
    D = D_tensor(voigt);
    D(4:6) = D(4:6) * 2;
    %}
    %
    % pass strain, stress and hardening variable to material subroutine
    % 3d voigt notation
    SIGMA_n = obj.stress_n(:, numgp);
    [C, SIGMA_np1, obj.hardening_np1(:,gp)] = mat.calcStressTangent(D, SIGMA_n, obj.hardening_n(:,gp));
    obj.stress(:, gp) = SIGMA_np1;
    % obj.strain(:, gp) = 0.5 * ( transpose(F_np1) * F_np1 - eye(dim) );
    %
    % compute rotation matrix
    % such that R'*eps*R = T*eps, R*sigma*R' = T'*sigma
    R_np1 = polarDecomp(F_np1);
    RT = transpose(R_np1);
    R1 = RT(:, [2,3,1]);
    R2 = RT([2,3,1], :);
    T22 = zeros(3, 3);
    for a = 1:3
        for b = 1:3
            temp = [1,2,3,1];
            c = temp(a+1);
            d = temp(b+1);
            T22(a, b) = RT(a, b) * RT(c, d) + RT(a, d) * RT(c, b);
        end
    end
    T = [RT.*RT, RT.*R1; 2*RT.*R2, T22];
    T(:, 4:6) = T(:, [5, 6, 4]);
    T(4:6, :) = T([5, 6, 4], :);
    %
    % compute material stiffness in current configuration
    sigma = transpose(T) * SIGMA_np1;
    Q = [2*sigma(1), 0,          0,          0,                    sigma(5),             sigma(6);
         0,          2*sigma(2), 0,          sigma(4),             0,                    sigma(6);
         0,          0,          2*sigma(3), sigma(4),             sigma(5),             0;
         0,          sigma(4),   sigma(4),   0.5*(sigma(2) + (3)), 0.5*sigma(6),         0.5*sigma(5);
         sigma(5),   0,          sigma(5),   0.5*sigma(6),         0.5*(sigma(1) + (3)), 0.5*sigma(4);
         sigma(6),   sigma(6),   0,          0.5*sigma(5),         0.5*sigma(4),         0.5*(sigma(1) + (2))];
    E = transpose(T) * C * T - Q;
%     E = C;
    E_2d = E([1,2,6], [1,2,6]);
    %
    % construct a 4*8 B matrix for finite deformation
    B = zeros(4, 8);
    B(1:2, 1:2:7) = dNdx;
    B(3:4, 2:2:8) = dNdx;
    %
    % form geometric stiffness term
    voigt32 = [1, 6; 6, 2];
    sigma_tensor = sigma(voigt32);
    sigma_matrix = zeros(4, 4);
    sigma_matrix(1:2, 1:2) = sigma_tensor;
    sigma_matrix(3:4, 3:4) = sigma_tensor;
    %
    % form material stiffness term
    temp = [1, 0, 0, 0; 0, 0, 0, 1; 0, 1, 1, 0];
    C_matrix = transpose(temp) * E_2d * temp;
    %
    % update Fint and K
    voigt34 = [1; 6; 6; 2];
    K = K + transpose(B) * (C_matrix + sigma_matrix) * B * det(Jcurr);
%     K = K + transpose(B) * (C_matrix) * B * det(Jref);
    Fint = Fint + transpose(B) * sigma(voigt34) * det(Jcurr);
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