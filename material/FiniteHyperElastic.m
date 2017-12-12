% finite strain hyperelastic material
% midterm exam in NLFEA class
% W(C) = 0.5 * miu * ( C - 3 ) - miu * log( J ) + 0.5 * lamda * ( J - 1 )^2
% Input: b = F^T * F (Voigt notation)
% Output: (1) Cauchy stress; (2) Material tangent (Voigt notation)
classdef FiniteHyperElastic < MaterialModel
    properties
        lamda
        miu
    end
    methods
        function obj = FiniteHyperElastic(inp)
            if nargin == 0
                obj.lamda = 40;
                obj.miu = 40;
            elseif nargin == 1
                if ~isequal(size(inp),[2,1])
                    error('Material parameters of Hyperelastic material shall be a 2x1 vector');
                end
                obj.lamda=inp(1);
                obj.miu=inp(2);
            else
                error('Invalid input for FiniteHyperElastic.m');
            end
        end
        function Initialize(obj, ~)
            obj.nhardening = 0; % no hardening variables for this model
        end
        function [D, sigma] = calcStressTangent(obj, eps)
            if isequal(size(eps),[3,1])
                J = sqrt( eps(1) * eps(2) - 0.25 * eps(3)^2 );
                if ~isreal(J)
                    error('Invalid strain for FiniteHyperElastic.m');
                end
                temp = obj.miu / J;
                D = [2*temp + obj.lamda, (2*J-1)*obj.lamda, 0;
                     (2*J-1)*obj.lamda, 2*temp + obj.lamda, 0;
                     0, 0, temp - (J-1)*obj.lamda];
                sigma = [temp * (eps(1) - 1) + (J-1)*obj.lamda;
                         temp * (eps(2) - 1) + (J-1)*obj.lamda;
                         0.5 * temp * eps(3)];
            elseif isequal(size(eps), [6,1])
                J = eps(1)*eps(2)*eps(3) + 0.25*eps(4)*eps(5)*eps(6);
                for temp = 1:3
                    J = J - 0.25*eps(temp)*eps(temp+3)^2;
                end
                J = sqrt(J);
                if ~isreal(J)
                    error('Invalid strain for FiniteHyperElastic.m');
                end
                temp = obj.miu / J;
                D = [2*temp + obj.lamda, (2*J-1)*obj.lamda, (2*J-1)*obj.lamda, 0, 0, 0;
                     (2*J-1)*obj.lamda, 2*temp + obj.lamda, (2*J-1)*obj.lamda, 0, 0, 0;
                     (2*J-1)*obj.lamda, (2*J-1)*obj.lamda, 2*temp + obj.lamda, 0, 0, 0;
                     0, 0, 0, temp - (J-1)*obj.lamda, 0, 0;
                     0, 0, 0, 0, temp - (J-1)*obj.lamda, 0;
                     0, 0, 0, 0, 0, temp - (J-1)*obj.lamda];
                sigma = [temp * (eps(1) - 1) + (J-1)*obj.lamda;
                         temp * (eps(2) - 1) + (J-1)*obj.lamda;
                         temp * (eps(3) - 1) + (J-1)*obj.lamda;
                         0.5 * temp * eps(4);
                         0.5 * temp * eps(5);
                         0.5 * temp * eps(6)];
            else
                error('strain must be Voigt notation')
            end
        end
    end
end