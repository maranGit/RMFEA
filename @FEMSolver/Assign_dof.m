function Assign_dof(obj, d)
for temp1 = 1:obj.nNodes
    for temp2 = 1:obj.ndofpn
        currDof = obj.nodes(temp1).ndof(temp2);
        if ~currDof.p
            currDof.v = d(currDof.pos);
        end
    end
end
end