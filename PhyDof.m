classdef PhyDof < handle
    properties
        p = 0 % (boolean) whether the dof is prescribed
        pos = 1 % (int) position in the global system (for free and prescribed)
        v = 0.0 % (double) value of dof, i.e. global displacement
        f = 0.0 % (double) force corresponding to dof
    end
    methods
        function obj = PhyDof(p, pos, v, f)
            if nargin == 0
                obj.p = 0;
                obj.pos = 1;
                obj.v = 0.0;
                obj.f = 0.0;
            elseif nargin == 4
                obj.p = p;
                obj.pos = pos;
                obj.v = v;
                obj.f = f;
            else
                error('Invalid initiation of PhyDof');
            end
        end
        function vUpdate(obj, deltaV)
            obj.v = obj.v + deltaV;
        end
        function fUpdate(obj, deltaF)
            obj.f = obj.f + deltaF;
        end
    end
end