classdef Hyper3d < PhyElement
% Ran Ma
% 12/8/2017
% 2D element subroutine for hyper-elastic-plastic material

    properties
        
    end
    %     properties (Dependent)
    %         J
    %     end
    methods
        function obj = Hyper3d()
            obj = obj@PhyElement();
        end
        function Initialize(obj, nhardening)
            % initialize hardening variables
            obj.n_hardening = nhardening;
            if nhardening == 0
                obj.hardening_n = 0;
                obj.hardening_np1 = 0;
            else
                obj.hardening_n = zeros(nhardening, 4);
                obj.hardening_np1 = zeros(nhardening, 4);
            end
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
            % updated Lagrange element subroutine for finite strain
            % update element stiffness and Fint for current displacement
            % Input: material subroutine, X, displacement
            % Output: (1) Fint; (2) element stiffness
            %
            % parameters for shape function
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
            disp = [obj.nedof.v]'; % should change to edofs in future version @@@
            % x1, x2, x3, x4
            % y2, y3, y4, y4
            x = X + reshape(disp, 3, []);
            %
            % loop over integration point for Fint and K
            Fint = zeros(obj.numDofs, 1);
            K = zeros(obj.numDofs, obj.numDofs);
            for gp = 1:lint
                %
                % compute gradient using FEAP shape function
                if nel == 4 || nel == 10 % tetrahedral element
                  [Wgt,ss] =  int3d_t(gp,lint,ib);
                  [~,shld,shls,be] = shltt(ss,nel,nel,der,bf);
                  [dNdX, ~, ~] = shgtt(X,nel,shld,shls,nel,bf,der,be);
                  [dNdx, ~, Jcurr] = shgtt(x,nel,shld,shls,nel,bf,der,be);
                else % brick element
                  [Wgt,ss] =  intpntb(gp,lint,ib);
                  [~,shld,shls,be] = shlb(ss,nel,nel,der,bf);
                  [dNdX, ~, ~] = shgb(X,nel,shld,shls,nel,bf,der,be);
                  [dNdx, ~, Jcurr] = shgb(x,nel,shld,shls,nel,bf,der,be);
                end
                %
                % deformation gradient
                F = x * dNdX; % deformation gradient
                b = F * transpose(F); % left Cauchy-Green strain tensor
                b_voigt = [b(1); b(5); b(9); 2*b(8); 2*b(7); 2*b(4)];
                [C, sigma] = mat.calcStressTangent(b_voigt);
                obj.stress(:, gp) = sigma;
                obj.strain(:, gp) = 0.5 * ( b_voigt - [1; 1; 1; 0; 0; 0] );
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
                K = K + Wgt * transpose(B) * C * B * Jcurr + Wgt * transpose(B2) * M * B2 * Jcurr;
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