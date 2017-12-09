classdef PhyElement < handle
    properties
        id
        lint % number of integration point
        neNodes % element nodes ( int )
        eNodes % element node vector ( vector<int> )
        nedof % element dof ( int )
        numDofs % number of dofs
        dofMap % map from element to global dofs ( vector<int> )
        matID % type of material in global
        ke % element stiffness matrix
        fde % force vector from essential BC
        foe % surface load and body force (and thermal force in future)
        fee % foe - fde
        strain % strain at each integration point
        stress % stress of step n+1
        stress_n % stress of step n
        Fint   % F internal at each integration point
        n_hardening   % number of hardening variables
        hardening_n   % hardening variables from step n ( vector<double> )
        hardening_np1 % hardening variables from step n+1 ( vector<double> )
    end
    methods
        function obj = PhyElement()
            % object constructor
            
        end
        function Update(obj)
            obj.stress_n = obj.stress;
            obj.hardening_n = obj.hardening_np1;
        end
        function output(obj, fid)
            sizeStress = size([obj.stress], 1);
            sizeStrain = size([obj.strain], 1);
            formatStress = repmat('%5.6e, ', 1, sizeStress);
            formatStress([end-1, end]) = '\n';
            formatStress = strcat('  ', formatStress);
            formatStrain = repmat('%5.6e, ', 1, sizeStrain);
            formatStrain([end-1, end]) = '\n';
            formatStrain = strcat('  ', formatStrain);
            fprintf(fid, '  stress at each integration point\n');
            fprintf(fid, formatStress, obj.stress);
            fprintf(fid, '  strain at each integration point\n');
            fprintf(fid, formatStrain, obj.strain);
%             fprintf(fid, '  consistent stiffness matrix at each integration point\n');
%             fprintf(fid, '%5.6e, %5.6e, %5.6e\n', obj.C);
        end
    end
    methods (Abstract)
        % Initializing hardening variables vector
        Initialize(obj, n_hardening)
        
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