classdef Small3d < PhyElement
% Ran Ma
% 12/10/2017
% small strain 3d element

    properties
        
    end
    methods
        function obj = Small3d()
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
            B = zeros(6, 3*nel); % B matrix
            for gp = 1:obj.lint
                % form shape function and global derivative
                X = reshape([obj.eNodes.coordinates], 3, []);
                if nel == 4 || nel == 10 % tetrahedral element
                  [Wgt,ss] =  int3d_t(gp,lint,ib);
                  [~,shld,shls,be] = shltt(ss,nel,nel,der,bf);
                  [dNdX, ~, Jdet] = shgtt(X,nel,shld,shls,nel,bf,der,be);
                else % brick element
                  [Wgt,ss] =  intpntb(gp,lint,ib);
                  [~,shld,shls,be] = shlb(ss,nel,nel,der,bf);
                  [dNdX, ~, Jdet] = shgb(X,nel,shld,shls,nel,bf,der,be);
                end
                %
                % form B matrix
                temp = zeros(6, 3);
                temp([1, 12, 17]) = 1;
                B(:, 1:3:end) = temp * transpose(dNdX);
                temp = zeros(6, 3);
                temp([6, 8, 16]) = 1;
                B(:, 2:3:end) = temp * transpose(dNdX);
                temp = zeros(6, 3);
                temp([5, 10, 15]) = 1;
                B(:, 3:3:end) = temp * transpose(dNdX);
                %
                % material subroutine
                obj.strain(:, gp) = B*disp;
                if obj.n_hardening == 0 % elastic model
                    [D, sigma] = mat.calcStressTangent(obj.strain(:, gp));
                else % plastic model
                    [D, sigma, obj.hardening_np1(:,gp)] = mat.calcStressTangent(obj.strain(:, gp), obj.hardening_n(:,gp));
                end
                obj.stress(:, gp) = sigma;
                fint = fint + Wgt * transpose(B) * sigma * Jdet;
                K = K + Wgt * transpose(B) * D * B * Jdet;
            end
            obj.ke = K;
            obj.Fint = fint;
            %
            % construct element force vector (including body force, surface loading and displacement force)
            % pDisp = [obj.nedof.v]; % collecte displacement of all dof
            % pDisp(obj.dofMap>0) = 0; % set displacement of free dof to 0
%             obj.fde = zeros(8, 1);
%             obj.fee = obj.foe - obj.fde;
        end
        
        Calculate_Stress_Strain(obj)
        
        SpecificOutput(obj)
        
    end
    
    % use static function to deal with shape function and numerical
    % integration
    methods (Static = true)
        
    end
end