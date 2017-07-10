function Calculate_ElementStiffness_Force(obj)
% calculate element stiffness (ke) and force vector (fde, foe and fee)
for e = 1:obj.ne %loop over elements
    obj.elements(e).Calculate_ElementStiffness_Force;
end
end