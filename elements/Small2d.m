classdef Small2d < PhyElement
% Ran Ma
% 12/8/2017
% small strain 2d element

    properties
        thickness = 1.0
        Bmat_1
        Bmat_2
        Bmat_3
        Bmat_4
    end
%     properties (Dependent)
%         J
%     end
    methods
        function obj = Small2d()
            obj = obj@PhyElement();
            obj.neNodes = 4;
            obj.nedof = 8;
            obj.fee = zeros(8, 1);
            obj.fde = zeros(8, 1);
            obj.foe = zeros(8, 1);
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
            %             obj.C = zeros(3, 1);
        end
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
            fint = zeros(obj.numDofs, 1);
            K = zeros(obj.numDofs, obj.numDofs);
            for gp = 1:obj.lint
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
                B = [dNdX(1, 1),          0, dNdX(2, 1),          0, dNdX(3, 1),          0, dNdX(4, 1),          0;
                              0, dNdX(1, 2),          0, dNdX(2, 2),          0, dNdX(3, 2),          0, dNdX(4, 2);
                     dNdX(1, 2), dNdX(1, 1), dNdX(2, 2), dNdX(2, 1), dNdX(3, 2), dNdX(3, 1), dNdX(4, 2), dNdX(4, 1)];
                 obj.strain(:, gp) = B*disp;
                 if obj.n_hardening == 0 % elastic model
                     [D, sigma] = mat.calcStressTangent(obj.strain(:,gp));
                 else % plastic model
                     [D, sigma, obj.hardening_np1(:,gp)] = mat.calcStressTangent(obj.strain(:,gp), obj.hardening_n(:,gp));
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