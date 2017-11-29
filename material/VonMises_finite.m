% Ran Ma
% 11/22/2017
%
% finite strain von Mises model with linear hardening
% project for NLFEA class
% based on Box 8.1 in <<computational inelasticity>>
% Input: F_n, F_np1, tau, hardening_n
% structure of old hardening variable: 
% 3d
%   1: alpha_n
%   2-7: beta_n
% Output: (1) stress; (2) tangent stiffness; (3) new hardening variables
classdef VonMises_finite < MaterialModel
    properties
        E
        nu
        K0
        K1
        H1
        alpha = 0.5 % mid-point rule
    end
    methods
        function obj = VonMises_finite(inp)
            if nargin == 0
                obj.E = 1200;
                obj.nu = 0.25;
                obj.K0 = 22;
                obj.K1 = 2.5;
                obj.H1 = 3.5;
            elseif nargin == 1
                if ~isequal(size(inp),[5,1])
                    error('Material parameters of Hyperelastic material shall be a 5x1 vector');
                end
                obj.E = inp(1);
                obj.nu = inp(2);
                obj.K0 = inp(3);
                obj.K1 = inp(4);
                obj.H1 = inp(5);
            else
                error('Invalid input for Linear Hyperelastic material');
            end
        end
        function Initialize(obj, dimension)
            if dimension == 2          % 2d von-Mises model
                obj.dim = 2;
                obj.nhardening = 7;
            elseif dimension == 3      % 3d von-Mises model
                obj.dim = 3;
                obj.nhardening = 7;
            else
                error('Invalid dimension for material subroutine');
            end
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
            if ~isequal(size(hardening_n), [7,1])
                error('hardening variable should be a 7x1 vector');
            end
            %
            % update SIGMA_np1 and hardening_np1
            [SIGMA_np1, hardening_np1] = obj.calcStress(D, SIGMA_n, hardening_n);
            %
            % compute consistent stiffness matrix
            % complex step derivative
            h = 1e-12;
            C = zeros(6,6);
            for k = 1:6
                A = D;
                if abs(D(k)) < h
                    hi = h;
                else
                    hi = h*abs(D(k));
                end
                A(k) = A(k) + 1i*hi;
                [SIGMA_inc, ~] = obj.calcStress(A, SIGMA_n, hardening_n);
                C(1:6,k) = 1/hi*imag(SIGMA_inc - SIGMA_np1);
            end
        end
        %
        % input complex D, output complex SIGMA_n
        function [SIGMA_np1, hardening_np1] = calcStress(obj, D, SIGMA_n, hardening_n)
            %
            % frequently used constants
            Diag = transpose(1:3);
            offDiag = transpose(4:6);
            I4 = diag([1.0, 1.0, 1.0, 0.5, 0.5, 0.5]); % fourth order identity tensor
            I2_tensor_I2 = zeros(6, 6); % delta tensor-products delta
            I2_tensor_I2(1:3, 1:3) = 1;
            %
            % Initialize material constants
            K_0 = obj.K0;
            K_1 = obj.K1;
            H_1 = obj.H1;
            mu = 0.5 * obj.E / (1 + obj.nu);
            lamda = obj.E * obj.nu / (1 + obj.nu) / (1 - 2 * obj.nu);
            kappa = lamda + 2 * mu / 3;
            %
            % Initialization
            ALPHA_n = hardening_n(1, 1);
            Q_n = hardening_n(2:7, 1);
            %
            % elastic predictor (isotropic)
            C_np1 = kappa * I2_tensor_I2 + 2 * mu * (I4 - I2_tensor_I2 / 3);
            SIGMA_tr = SIGMA_n + C_np1 * D;
            Q_tr = Q_n;
            ksi_tr = SIGMA_tr - Q_tr;
            ksi_tr(Diag) = ksi_tr(Diag) - sum(ksi_tr(Diag)) * ones(3, 1) / 3;
            temp1 = ksi_tr .* conj(ksi_tr);
            temp2 = ksi_tr(offDiag) .* conj(ksi_tr(offDiag));
            normKsiTr = sqrt(sum(temp1) + sum(temp2));
            f_tr = normKsiTr - sqrt(2/3) * (K_0 + K_1 * ALPHA_n);
            %
            % radial return method
            if f_tr <= 0
                SIGMA_np1 = SIGMA_tr;
                hardening_np1 = hardening_n;
            else
                dgamma = f_tr / (2 * mu + 2 * (H_1 + K_1) / 3);
                N_np1 = ksi_tr / normKsiTr;
                SIGMA_np1 = SIGMA_tr - 2 * mu * dgamma * N_np1;
                Q_np1 = Q_tr + 2 * dgamma * H_1 * N_np1 / 3;
                ALPHA_np1 = ALPHA_n + sqrt(2/3) * dgamma;
                hardening_np1 = [ALPHA_np1; Q_np1];
            end
            
        end
    end
end