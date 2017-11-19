% general material model class
% Input: epsilon_ij, properties, data(plastic strain, flow stress, ...)
% Output: sigma_ij, C_ijkl
classdef MaterialModel < handle
    properties
        sigma % stress, Voigt notation
        C     % tangent stiffness, Voigt notation
        dim   % dimension of the model, normally either 2 or 3
        nhardening  % number of internal hardening variables
    end
    methods
        function obj = MaterialModel()
            
        end
    end
    methods (Abstract)
        calcStressTangent(obj)
        Initialize(obj)
    end
end