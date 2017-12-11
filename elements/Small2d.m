classdef Small2d < PhyElement
% Ran Ma
% 12/8/2017
% small strain 2d element

    properties
        thickness = 1.0
    end
    methods
        function obj = Small2d()
            obj = obj@PhyElement();
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
                obj.hardening_n = zeros(nhardening, obj.lint);
                obj.hardening_np1 = zeros(nhardening, obj.lint);
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
            % update element stiffness and Fint for current displacement
            % Input: nodal displacement
            % Output: (1) residual; (2) element stiffness
            %
            nel = obj.neNodes;
            lint = obj.lint;
            der = 0;
            bf = 0;
            ib = 0;
            disp = [obj.nedof.v]'; % should change to edofs in future version @@@
            fint = zeros(obj.numDofs, 1); % internal force
            K = zeros(obj.numDofs, obj.numDofs); % consistent stiffness
            B = zeros(3, 2*nel); % B matrix
            for gp = 1:obj.lint
                % form shape function and global derivative
                X = reshape([obj.eNodes.coordinates], 2, []);
                if nel == 3 || nel == 6
                    [Wgt,litr,lits] =  intpntt(gp,lint,ib);
                    [~,shld,shls,be] = shlt(litr,lits,nel,nel,der,bf);
                    [dNdX, ~, Jdet] = shgt(X,nel,shld,shls,nel,bf,der,be);
                else
                    [Wgt,litr,lits] =  intpntq(gp,lint,ib);
                    [~,shld,shls,be] = shlq(litr,lits,nel,nel,der,bf);
                    [dNdX, ~, Jdet] = shgq(X,nel,shld,shls,nel,bf,der,be);
                end
                %
                % form B matrix
                B(:, 1:2:end) = [1, 0; 0, 0; 0, 1] * transpose(dNdX);
                B(:, 2:2:end) = [0, 0; 0, 1; 1, 0] * transpose(dNdX);
                %
                % material subroutine
                obj.strain(:, gp) = B*disp;
                epsilon = zeros(6, 1);
                epsilon([1,2,6]) = obj.strain(:, gp);
                if obj.n_hardening == 0 % elastic model
                    [D, sigma] = mat.calcStressTangent(epsilon);
                else % plastic model
                    [D, sigma, obj.hardening_np1(:,gp)] = mat.calcStressTangent(epsilon, obj.hardening_n(:,gp));
                end
                obj.stress(:, gp) = sigma;
                fint = fint + Wgt * transpose(B) * sigma([1,2,6]) * Jdet;
                K = K + Wgt * transpose(B) * D([1,2,6],[1,2,6]) * B * Jdet;
            end
            obj.ke = K;
            obj.Fint = fint;
            %
            % construct element force vector (including body force, surface loading and displacement force)
            % pDisp = [obj.nedof.v]; % collecte displacement of all dof
            % pDisp(obj.dofMap>0) = 0; % set displacement of free dof to 0
            obj.fde = zeros(8, 1);
            obj.fee = obj.foe - obj.fde;
        end
        
        Calculate_Stress_Strain(obj)
        
        SpecificOutput(obj)
        
    end
    
    % use static function to deal with shape function and numerical
    % integration
    methods (Static = true)
        
    end
end