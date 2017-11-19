% small strain hyperelastic material
% HW2 in NLFEA class
% W(eps) = a*eps_mm*log(1+eps_nn)/2 + 3*b/*eps_mn*eps_mn/2
% Input: strain (Voigt notation)
% Output: (1) stress; (2) tangent stiffness. (Voigt notation)
classdef Hyperelastic < MaterialModel
    properties
        a
        b
    end
    methods
        function obj = Hyperelastic(inp)
            if nargin == 0
                obj.a = 40;
                obj.b = 60;
            elseif nargin == 1
                if ~isequal(size(inp),[2,1])
                    error('Material parameters of Hyperelastic material shall be a 2x1 vector');
                end
                obj.a=inp(1);
                obj.b=inp(2);
            else
                error('Invalid input for Linear Hyperelastic material');
            end
        end
        function Initialize(obj, ~)
            obj.nhardening = 0; % no hardening variables for this model
        end
        function [D, sigma] = calcStressTangent(obj, eps)
            if isequal(size(eps),[3,1])
                traceEps = eps(1) + eps(2);
                temp = 0.5*obj.a*(2+traceEps)/(1+traceEps)^2;
                D = [temp+3*obj.b, temp, 0;
                     temp, temp+3*obj.b, 0;
                     0, 0, 1.5*obj.b];
                temp = 0.5 * obj.a * (log(1+traceEps) + traceEps/(1+traceEps));
                sigma = [temp + 3*obj.b*eps(1);
                         temp + 3*obj.b*eps(2);
                         1.5 * obj.b * eps(3)];
            elseif isequal(size(eps), [6,1])
                D = 0;
                sigma = 0;
            else
                error('strain must be Voigt notation')
            end
        end
    end
end