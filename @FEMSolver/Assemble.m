function Assemble(obj)
% assemble stiffness matrix and essential BC caused force
% step 1: initialize K
obj.K = zeros(obj.nf, obj.nf);
for e = 1:obj.ne %loop over elements
    % fee = feo; %element total force = element all forces except essential force
    currEle = obj.elements(e);
    nedof = size(currEle.nedof, 1);
    dofMap = currEle.dofMap;
    for i = 1:nedof %loop over rows of ke; nedof = element # dof
        I = dofMap(i); %local to global dof map Me t
        if (I > 0) %I corresponds to a free dof, we skip prescribed dofs
            for j = 1:nedof %loop over columns of ke
                J = dofMap(j); %global dof corresponding to j
                if (J > 0) %now both I and J are free and can add ke(i,j) to global K
                    obj.K(I, J) = obj.K(I, J) + currEle.ke(i, j);
                else %J < 0, prescribed dof j; add contributions of fDe = keae to fee
                    % currEle.fee(i) = currEle.fee(i) - ke(i, j) * edofs(j); %edofs: element dofs = ae
                end
            end
            obj.F(I) = obj.F(I) + currEle.fee(i); % element’s total force fee component i’th is computed->added to F(I)
        end
    end
end
end