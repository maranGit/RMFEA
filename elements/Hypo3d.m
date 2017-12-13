classdef Hypo3d < PhyElement
    % Ran Ma
    % 11/22/2017
    % generalized element subroutine for all finite strain plastic model
    % element subroutine is based on warp3d
    
    properties
        thickness = 1.0
        tau_n % kirchoff stress of the previous step
        alpha = 0.5; % alpha configuration for mid-point integration
        stress_np1
    end
    methods
        function obj = Hypo3d()
            obj = obj@PhyElement();
            obj.stress = zeros(6, obj.lint);
            obj.strain = zeros(6, obj.lint);
            %obj.eNodes(4, 1) = PhyNode();
            %obj.edofs(8, 1) = PhyDof();
        end
        function Initialize(obj, nhardening)
            % initialize hardening variables
            obj.n_hardening = nhardening;
            if nhardening == 0
                obj.hardening_n = zeros(1, obj.lint);
                obj.hardening_np1 = zeros(1, obj.lint);
            else
                obj.hardening_n = zeros(nhardening, obj.lint);
                obj.hardening_np1 = zeros(nhardening, obj.lint);
            end
            obj.stress_n = zeros(6, obj.lint);
            obj.stress = zeros(6, obj.lint);
            obj.fee = zeros(obj.numDofs, 1);
            obj.fde = zeros(obj.numDofs, 1);
            obj.foe = zeros(obj.numDofs, 1);
            obj.ke = zeros(obj.numDofs, obj.numDofs);
            obj.Fint = zeros(obj.numDofs, 1);
        end
        % detail of virtual function in PhyElement.m
        % get surface area of the element
        setGeometry(obj)
        
        % read material properties, and change to plane strain
        setInternalMaterialProperties(obj)
        
        % form element stiffness matrix and force vector
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
            % parameters for FEAP shape function
            nel = obj.neNodes;
            lint = obj.lint;
            der = 0;
            bf = 0;
            ib = 0;
            %
            % X1, X2, X3, X4
            % Y2, Y3, Y4, Y4
            X = reshape([obj.eNodes.coordinates], 3, []);
            % u1, v1, u2, v2, u3, v3, u4, v4
            disp_n = [obj.nedof.vn]'; % should change to edofs in future version @@@
            disp_np1 = [obj.nedof.v]';
            numDof = length(disp_n);
            % x1, x2, x3, x4
            % y2, y3, y4, y4
            x_n = X + reshape(disp_n, 3, []);
            x_np1 = X + reshape(disp_np1, 3, []);
            %
            % loop over integration point for Fint and K
            Fint = zeros(numDof, 1);
            K = zeros(numDof, numDof);
            for gp = 1:lint
                % compute gradient using FEAP shape function
                if nel == 4 || nel == 10 % tetrahedral element
                  [Wgt,ss] =  int3d_t(gp,lint,ib);
                  [~,shld,shls,be] = shltt(ss,nel,nel,der,bf);
                  [dNdX, ~, ~] = shgtt(X,nel,shld,shls,nel,bf,der,be);
                  [dNdx, ~, Jcurr] = shgtt(x_np1,nel,shld,shls,nel,bf,der,be);
                else % brick element
                  [Wgt,ss] =  intpntb(gp,lint,ib);
                  [~,shld,shls,be] = shlb(ss,nel,nel,der,bf);
                  [dNdX, ~, ~] = shgb(X,nel,shld,shls,nel,bf,der,be);
                  [dNdx, ~, Jcurr] = shgb(x_np1,nel,shld,shls,nel,bf,der,be);
                end
                %
                % compute F_np1
                F_n = x_n * dNdX; % deformation gradient
                F_np1 = x_np1 * dNdX;
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
                SIGMA_n = obj.stress_n(:, gp);
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
                %
                % form B matrix for material stiffness
                B = zeros(6, 3*nel);
                temp = zeros(6, 3);
                temp([1, 12, 17]) = 1;
                B(:, 1:3:end) = temp * transpose(dNdx);
                temp = zeros(6, 3);
                temp([6, 8, 16]) = 1;
                B(:, 2:3:end) = temp * transpose(dNdx);
                temp = zeros(6, 3);
                temp([5, 10, 15]) = 1;
                B(:, 3:3:end) = temp * transpose(dNdx);
                %
                % form B2 matrix for geometric stiffness
                B2 = zeros(9, 3*nel);
                temp = zeros(9, 3);
                temp([1, 13, 25]) = 1;
                B2(:, 1:3:end) = temp * transpose(dNdx);
                temp = zeros(9, 3);
                temp([2, 14, 26]) = 1;
                B2(:, 2:3:end) = temp * transpose(dNdx);
                temp = zeros(9, 3);
                temp([3, 15, 27]) = 1;
                B2(:, 3:3:end) = temp * transpose(dNdx);
                %
                % form geometric stiffness term
                M = [sigma(1)*eye(3), sigma(6)*eye(3), sigma(5)*eye(3);
                    sigma(6)*eye(3), sigma(2)*eye(3), sigma(4)*eye(3);
                    sigma(5)*eye(3), sigma(4)*eye(3), sigma(3)*eye(3)];
                %
                % update Fint and K
                K = K + Wgt * transpose(B) * E * B * Jcurr + Wgt * transpose(B2) * M * B2 * Jcurr;
                Fint = Wgt * Fint + transpose(B) * sigma * Jcurr;
            end
            obj.ke = K;
            obj.Fint = Fint;
            %
            % construct element force vector (including body force, surface loading and displacement force)
            % pDisp = transpose(disp); % collect displacement of all dof
            % pDisp(obj.dofMap>0) = 0; % set displacement of free dof to 0
            % obj.fde = K * transpose(pDisp);
            obj.fde = zeros(obj.numDofs, 1);
            obj.fee = obj.foe - obj.fde;
            %
        end
        
        Calculate_Stress_Strain(obj)
        
        SpecificOutput(obj)
        
    end
    
    % use static function to deal with shape function and numerical
    % integration
    methods (Static = true)
        
    end
end