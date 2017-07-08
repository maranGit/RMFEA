function setPositions_F(obj)
Force = zeros(obj.nf, 1);
temp = 1;
for temp2 = 1:obj.nf
    Force(temp2, 1) = obj.dofs(temp2).f;
    temp = temp + 1;
end
obj.F = Force;
end