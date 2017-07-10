classdef PhyElement < handle
    properties
        id
        neNodes % element nodes ( int )
        eNodes % element node vector ( vector<int> )
        nedof % element dof ( int )
        edofs % vector of dofs
        dofMap % map from element to global dofs ( vector<int> )
        matID % type of material in global
        ke % element stiffness matrix
        fde % force vector from essential BC
        foe % surface load and body force (and thermal force in future)
        fee % foe - fde
        strain % strain at each integration point
        stress % stress at each integration point
    end
    methods
        function obj = PhyElement()
            % object constructor
            
        end
    end
    methods (Abstract)
        % coordinate -> geometry
        setGeometry(obj)
        
        % get material properties from global
        setInternalMaterialProperties(obj)
        
        % compute element stiffness/force
        Calculate_ElementStiffness_Force(obj)
        % 1. calculate foe: sum of all forces, but fde
        % 2. calculate fde = ke * ae
        % 3. calculate fee = fde + foe
        
        Calculate_Stress_Strain(obj)
        
        SpecificOutput(obj)
        
    end
end