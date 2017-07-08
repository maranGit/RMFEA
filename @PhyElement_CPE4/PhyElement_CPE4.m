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
        E = 30.0e6
        nu = 0.3
    end
%     properties (Dependent)
%         J
%     end
    methods
        function obj = PhyElement_CPE4()
            obj.neNodes = 4;
            obj.nedof = 8;
            obj.fee = zeros(8, 1);
            %obj.eNodes(4, 1) = PhyNode();
            %obj.edofs(8, 1) = PhyDof();
        end
        J = formJ(obj)
        
        % detail of virtual function in PhyElement.m
        % get surface area of the element
        setGeometry(obj)
        
        % read material properties, and change to plane strain
        setInternalMaterialProperties(obj)
        
        % form element stiffness matrix and force vector
        Calculate_ElementStiffness_Force(obj)
        
        SpecificOutput(obj)
        
    end
    
    % use static function to deal with shape function and numerical
    % integration
    methods (Static = true)
        
    end
end