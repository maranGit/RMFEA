% small strain hyperelastic material
% HW2 in NLFEA class
% W(eps) = a*eps_mm*log(1+eps_nn)/2 + 3*b/*eps_mn*eps_mn/2
% Input: strain (Voigt notation)
% Output: (1) stress; (2) tangent stiffness. (Voigt notation)
classdef Hypoelastic < MaterialModel
    properties
        E
        nu
    end
    methods
        function obj = Hypoelastic(inp)
            if nargin == 0
                obj.E = 1200;
                obj.nu = 0.25;
            elseif nargin == 1
                if ~isequal(size(inp),[2,1])
                    error('Material parameters of Hypoelastic material shall be a 2x1 vector');
                end
                obj.E = inp(1);
                obj.nu = inp(2);
            else
                error('Invalid input for Linear Hypoelastic material');
            end
        end
        function Initialize(obj, ~)
            obj.nhardening = 0;
        end
        function [C, SIGMA_np1, hardening_np1] = calcStressTangent(obj, D, SIGMA_n, hardening_n)
            %
            % check input
            if ~isequal(size(D), [6,1])
                error('D should be a 6x1 strain vector');
            end
            if ~isequal(size(SIGMA_n), [6,1])
                error('SIGMA_n should be a 6x1 stress vector');
            end
            %
            % frequently used constants
            I4 = diag([1.0, 1.0, 1.0, 0.5, 0.5, 0.5]); % fourth order identity tensor
            I2_tensor_I2 = zeros(6, 6); % delta tensor-products delta
            I2_tensor_I2(1:3, 1:3) = 1;
            %
            % Initialize material constants
            mu = 0.5 * obj.E / (1 + obj.nu);
            lamda = obj.E * obj.nu / (1 + obj.nu) / (1 - 2 * obj.nu);
            kappa = lamda + 2 * mu / 3;
            %
            % elastic predictor (isotropic)
            C = kappa * I2_tensor_I2 + 2 * mu * (I4 - I2_tensor_I2 / 3);
            SIGMA_np1 = SIGMA_n + C * D;
            hardening_np1 = hardening_n;
        end
    end
end