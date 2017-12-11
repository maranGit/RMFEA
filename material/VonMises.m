% small strain von Mises model
% HW4 in NLFEA class
% beta_dot = 2/3 * gamma * H1 * eta / norm(eta)
% alpha_dot = sqrt(2/3) * gamma
% Input: strain (Voigt notation), old hardening variables
% structure of old hardening variable: 
% 2d
%   1-3: ep_n
%   4: alpha_n
%   5-7: beta_n
% 3d
%   1-6: ep_n
%   7: alpha_n
%   8-13: beta_n
% Output: (1) stress; (2) tangent stiffness; (3) new hardening variables
classdef VonMises < MaterialModel
    properties
        E
        nu
        K0
        K1
        H1
    end
    methods
        function obj = VonMises(inp)
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
                obj.nhardening = 13;
            elseif dimension == 3      % 3d von-Mises model
                obj.nhardening = 13;
            else
                error('Invalid dimension for material subroutine');
            end
        end
        function [C, sigma, hardening_np1] = calcStressTangent(obj, eps, hardening_n)
%             voigt2 = [1;5;4];
%             voigt2_inv = [1, 3; 3, 2];
            voigt3 = [1; 5; 9; 8; 7; 4];
            voigt3_inv = [1, 6, 5; 6, 2, 4; 5, 4, 3];
            I4 = diag([1.0, 1.0, 1.0, 0.5, 0.5, 0.5]); % fourth order identity tensor
            I2_tensor_I2 = zeros(6, 6); % delta tensor-products delta
            I2_tensor_I2(1:3, 1:3) = 1;
            off_diag = [2; 3; 4; 6; 7; 8]; % off-diagonal element of 3x3 matrix
            %
            % Initialization
            ep_n = hardening_n(voigt3_inv);
            ep_n(off_diag) = ep_n(off_diag) * 0.5;
            alpha_n = hardening_n(7, 1);
            beta_n = hardening_n(voigt3_inv + 7);
            %
            % strain tensor
            if ~isequal(size(eps), [6,1])
                error('strain must be Voigt notation (6x1) vector');
            end
%             dim = 3;
            eps_new = eps(voigt3_inv);
            eps_new(off_diag) = eps_new(off_diag) * 0.5;
            eps = eps_new;
            %
            % Initialize
            K_0 = obj.K0;
            K_1 = obj.K1;
            H_1 = obj.H1;
            mu = 0.5 * obj.E / (1 + obj.nu);
            lamda = obj.E * obj.nu / (1 + obj.nu) / (1 - 2 * obj.nu);
            kappa = lamda + 2 * mu / 3;
            %
            % compute trial elastic stress
            e_np1 = eps - (trace(eps) * eye(3))/3;
            s_tr = 2 * mu * (e_np1 - ep_n);
            ksi_tr = s_tr - beta_n;
            %
            % check yield condition
            normKsiTr = sqrt(sumsqr(ksi_tr));
            f_tr = normKsiTr - sqrt(2/3) * (K_0 + K_1 * alpha_n);
            if f_tr <= 0
                sigma = s_tr + kappa * trace(eps) * eye(3);
                C_np1 = kappa * I2_tensor_I2 + 2 * mu * (I4 - I2_tensor_I2 / 3);
                hardening_np1 = hardening_n;
                sigma = sigma(voigt3);
                C = C_np1;
                return
            end
            %
            % compute n_np1 and find dGamma
            [dGamma, alpha_np1] = obj.gamma(normKsiTr, alpha_n);
            n_np1 = ksi_tr / normKsiTr;
            %
            % update back stress, plastic strain and stress
            beta_np1 = beta_n + sqrt(2/3) * H_1 * (alpha_np1 - alpha_n) * n_np1;
            ep_np1 = ep_n + dGamma * n_np1;
            sigma_np1 = s_tr + kappa * trace(eps) * eye(3) - 2 * mu * dGamma * n_np1;
            %
            % compute consistent elastoplastic tangent moduli
            theta = 1 - 2 * mu * dGamma / normKsiTr;
            theta_bar = 1 / (1 + (K_1 + H_1) / 3 / mu) - (1 - theta);
            n_tensor_n = n_np1(voigt3) * transpose(n_np1(voigt3));
            C_np1 = kappa * I2_tensor_I2 + 2 * mu * theta * (I4 - I2_tensor_I2 / 3) - 2 * mu * theta_bar * n_tensor_n;
            %
            % transfer to Voigt notation
            sigma = sigma_np1(voigt3);
            C = C_np1;
            hardening_np1 = [ep_np1(voigt3); alpha_np1; beta_np1(voigt3)];
            hardening_np1(4:6) = hardening_np1(4:6) * 2;
        end
        function [dGamma, alpha_np1] = gamma(obj, ksi_tr, alpha_n)
            % Initialize
            tol = 1e-12;
            dGamma = 0;
            alpha_np1 = alpha_n;
            K_0 = obj.K0;
            K_1 = obj.K1;
            H_1 = obj.H1;
            mu = 0.5 * obj.E / (1 + obj.nu);
            % Iterate
            g1 = -sqrt(2/3) * (K_0 + K_1 * alpha_np1) + ksi_tr;
            g2 = 2*mu*dGamma + sqrt(2/3)*H_1*(alpha_np1-alpha_n);
            g = g1 - g2;
            while abs(g) > tol
                dg = -2 * mu * (1 + (H_1 + K_1) / 3 / mu);
                dGamma = dGamma - g / dg;
                g1 = -sqrt(2/3) * (K_0 + K_1 * alpha_np1) + ksi_tr;
                g2 = 2*mu*dGamma + sqrt(2/3)*H_1*(alpha_np1-alpha_n);
                g = g1 - g2;
                alpha_np1 = alpha_n + sqrt(2/3) * dGamma;
            end
        end
    end
end