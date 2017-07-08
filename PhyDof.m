classdef PhyDof < handle
    properties
        p = 0 % (boolean) whether the dof is prescribed
        pos = 1 % (int) position in the global system (for free and prescribed)
        v = 0.0 % (double) value of dof, i.e. global displacement
        f = 0.0 % (double) force corresponding to dof
    end
    methods
        
    end
end