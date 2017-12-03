classdef PhyElement_CPE4 < PhyElement
% plane strain element with full integration
%           ^ ksi 2
%           |
% 4---------+----------3
% |   *     |      *   |
% |         |          |
% |---------+---------------> ksi 1
% |         |          |
% |   *     |      *   |
% 1---------+----------2
%           |
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
        function obj = PhyElement_CPE4()
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
        J = formJ(obj)
        
        % detail of virtual function in PhyElement.m
        % get surface area of the element
        setGeometry(obj)
        
        % read material properties, and change to plane strain
        setInternalMaterialProperties(obj)
        
        % form element stiffness matrix and force vector
        Calculate_ElementStiffness_Force(obj, mat)
        
        Calculate_Stress_Strain(obj)
        
        SpecificOutput(obj)
        
    end
    
    % use static function to deal with shape function and numerical
    % integration
    methods (Static = true)
        
    end
end