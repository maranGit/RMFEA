classdef PhyNode < handle
    properties
        id
        nndof % (int) number of dofs
        ndof % vector <PhyDof>
        coordinates % vector
    end
    methods
        function obj = PhyNode(id, nndof, coordinates)
            if nargin == 0
                obj.id = 1;
                obj.nndof = 1;
                obj.coordinates = 0;
            elseif nargin == 3
                obj.id = id;
                obj.nndof = nndof;
                obj.coordinates = coordinates;
            else
                error('Invalid input in PhyNode')
            end
        end
        
        set_nndof(obj)
        
        UpdateNodePrescribedDofForces(obj)
        
    end
end