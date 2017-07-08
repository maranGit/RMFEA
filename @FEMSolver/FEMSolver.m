classdef FEMSolver < handle
    properties
        dim % int
        ndofpn % int
        nNodes % (int) number of nodes in the domain
        ne % number of elements ( int )
        nodes % vector<PhyNode>
        elements % vector<PhyElement*>
        dofs % vector of global free dofs
        nf % number of free dofs
        np % number of prescribed dofs
        ndof % ndof = nf + np
        K % stiffness matrix
        F % Force vector (external load)
        Fp % vector of prescribed dof forces
        nmats % number of materials
    end
    methods
        function obj = FEMSolver(dim)		% works as function with one or zero argument
            if nargin == 0
                obj.dim = 2;
                obj.ndofpn = 2;
            else
                obj.dim = dim;
                obj.ndofpn = dim;
            end
        end
        
        % All steps together:
        function FEMSolve(obj, runName)
            obj.Input(runName);
            %obj.Calculate_ElementStiffness_Force();
            %obj.Assemble();
        end
    end
    
    methods (Access = private)
        Input(obj, runName)
        
        % Step 3: compute nf from ndof and np, initialize stiffness matrix and force vector
        setSizes(obj)
        
        % Step 4: set prescribed dofs: already done when reading the input file
        % Step 5: Set global free nodal dof: already done when reading the input file
        % Step 6 and Step 7: dof positions; Step 7: Set F
        setPositions_F(obj)
        
        % Step 8: Element dof maps Me
        % Step 9: Set element dofs ae
        setElementDofMap_ae(obj)
        
        % Step 10: Compute element stiffness
        Calculate_ElementStiffness_Force(obj)
        
        % Step 11: Assembly from local to global system
        Assemble(obj)
        
        % Step 12: Solve global (free) dof a from Ka = F
        % successful solution returns true
        Solve_Dofs(obj)
        
        % Step 13: Assign a to nodes and elements
        Assign_dof(obj)
        
        % Step 14: Compute prescribed dof forces
        UpdateFpNodalPrescribedForces(obj)
    end
end