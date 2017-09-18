% general material model class
% Input: epsilon_ij, properties, data(plastic strain, flow stress, ...)
% Output: sigma_ij, C_ijkl
classdef MaterialModel < handle
    properties
        sigma % stress, Voigt notation
        C     % tangent stiffness, Voigt notation
    end
    methods
        function obj = MaterialModel()
                    
        end
    end
    methods (Abstract)
        calcStressTangent(obj)
    end
end