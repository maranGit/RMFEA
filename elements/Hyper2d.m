classdef Hyper2d < PhyElement
% Ran Ma
% 12/8/2017
% 2D element subroutine for hyper-elastic-plastic material

    properties
        thickness = 1.0
    end
    %     properties (Dependent)
    %         J
    %     end
    methods
        function obj = Hyper2d()
            obj = obj@PhyElement();
            obj.neNodes = 4;
            obj.nedof = 8;
            obj.fee = zeros(8, 1);
            obj.fde = zeros(8, 1);
            obj.foe = zeros(8, 1);
            obj.stress = zeros(3, 4);
            obj.strain = zeros(3, 4);
            %obj.eNodes(4, 1) = PhyNode();
            %obj.edofs(8, 1) = PhyDof();
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
            for gp = 1:lint
                %
                % compute gradient using FEAP shape function
                if nel == 3 || nel == 6
                    [Wgt,litr,lits] =  intpntt(gp,lint,ib);
                    [~,shld,shls,be] = shlt(litr,lits,nel,nel,der,bf);
                    [dNdX, ~, ~] = shgt(X,nel,shld,shls,nel,bf,der,be);
                    [dNdx, ~, Jcurr] = shgt(x,nel,shld,shls,nel,bf,der,be);
                else
                    [Wgt,litr,lits] =  intpntq(gp,lint,ib);
                    [~,shld,shls,be] = shlq(litr,lits,nel,nel,der,bf);
                    [dNdX, ~, ~] = shgq(X,nel,shld,shls,nel,bf,der,be);
                    [dNdx, ~, Jcurr] = shgq(x,nel,shld,shls,nel,bf,der,be);
                end
                %
                % deformation gradient
                F = x * dNdX; % deformation gradient
                % construct a 4*8 B matrix for finite deformation
                B = zeros(4, 2*nel);
                B(1:2, 1:2:end) = transpose(dNdx);
                B(3:4, 2:2:end) = transpose(dNdx);
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
                K = K + Wgt * transpose(B) * (C_matrix + sigma_matrix) * B * Jcurr;
                Fint = Wgt * Fint + transpose(B) * [sigma(1); sigma(3); sigma(3); sigma(2)] * Jcurr;
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
        
        Calculate_Stress_Strain(obj)
        
        SpecificOutput(obj)
        
    end
    
    % use static function to deal with shape function and numerical
    % integration
    methods (Static = true)
        
    end
end