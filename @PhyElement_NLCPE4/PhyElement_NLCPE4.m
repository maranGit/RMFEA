classdef PhyElement_NLCPE4 < PhyElement
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
        function obj = PhyElement_NLCPE4()
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