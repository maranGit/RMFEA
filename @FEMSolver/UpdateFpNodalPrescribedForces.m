function UpdateFpNodalPrescribedForces(obj)
% Ran Ma
% 2017/07/09
% Loop over all the elements to calculate nodel force
% and stress & strain at each integration point
% Will compute reaction force in the future
for temp = 1:obj.ne
    obj.elements(temp).Calculate_Stress_Strain;
end
end