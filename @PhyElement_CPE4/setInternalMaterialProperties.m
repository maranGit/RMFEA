function setInternalMaterialProperties(obj)
% D matrix
E = obj.E;
nu = obj.nu;

% plain strain
E1 = E * (1 - nu) / (1 + nu) / ( 1 - 2 * nu );
E2 = nu * E1 / (1 - nu);
G = 0.5 * E / ( 1 + nu );

% plain stress
% E1 = E / (1 - nu ^ 2);
% E2 = nu * E1;
% G = 0.5 * E / ( 1 + nu );

obj.E1 = E1;
obj.E2 = E2;
obj.G = G;
end