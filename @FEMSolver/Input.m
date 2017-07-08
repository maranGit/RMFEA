function Input(obj, runName)
%% read from input file
if exist('runName', 'var')
    run(runName);
else
    coordinates = [0,0; 1,0; 0,1; 1,1; 2,0; 2,1];
    incidences = [1 1 2 4 3; 2 2 5 6 4];
    eType = 'PhyElement_CPE4';
    EBC = [1, 1, 0;
        1, 2, 0;
        3, 1, 0];
    NBC = [5, 1, 1;
        6, 1, 1];
end

%% check input file
nNodes_temp = size(coordinates, 1); % int, total nodes
ne_temp = size(incidences, 1); % int, total elements
ndof_temp = nNodes_temp * obj.dim; % int, total dofs
np_temp = size(EBC, 1); % int, total prescribed dofs
nf_temp = ndof_temp - np_temp; % int, total free dofs

%% initialization
% initialize problem size
obj.nNodes = nNodes_temp;
obj.ne = ne_temp;
obj.ndof = ndof_temp;
obj.np = np_temp;
obj.nf = nf_temp;
dofAll(ndof_temp, 1) = PhyDof();
dofP = (EBC(:,1) - 1) * obj.dim + EBC(:,2); % (int), prescribed dof list
dofF = (NBC(:,1) - 1) * obj.dim + NBC(:,2); % (int), loaded dof list
dofAssigned = zeros(ndof_temp, 1); % which dof has essential BC

% initialize prescribed dofs
% sequence: 1x, 1y, 1z, 2x, 2y, 2z, ...
for temp = 1:np_temp
    currDof = dofP(temp);
    dofAll(currDof).p = 1;
    dofAll(currDof).pos = -temp;
    dofAll(currDof).v = EBC(temp, 3);
    dofAll(currDof).f = 0;
    dofAssigned(currDof) = 1;
end

% initialize free dofs
temp2 = 1;
for temp = 1:ndof_temp
    if dofAssigned(temp) == 0
        dofAll(temp).p = 0;
        dofAll(temp).pos = temp2;
        dofAll(temp).v = 0;
        dofAll(temp).f = 0;
        temp2 = temp2 + 1;
    end
end
for temp = 1:size(NBC, 1)
    currDof = dofF(temp);
    dofAll(currDof).p = 0;
    dofAll(currDof).v = 0;
    dofAll(currDof).f = NBC(temp, 3);
end
obj.dofs = dofAll(dofAssigned == 0);

% initialize nodes
curDim = obj.dim;
nodeAll(nNodes_temp, 1) = PhyNode();
for temp = 1:nNodes_temp
    nodeAll(temp).id = temp;
    nodeAll(temp).nndof = curDim;
    nodeAll(temp).ndof = dofAll((temp*curDim-curDim+1):(temp*curDim));
    nodeAll(temp).coordinates = coordinates(temp, :);
end
obj.nodes = nodeAll;

% initialize elements
eleAll(ne_temp, 1) = eval(eType); % avoid using 'switch'
for temp = 1:ne_temp
    eleAll(temp).id = incidences(temp, 1);
    eleAll(temp).eNodes = obj.nodes(incidences(temp, 2:end));
    eleAll(temp).nedof = reshape([eleAll(temp).eNodes.ndof], [], 1);
    eleAll(temp).dofMap = [eleAll(temp).nedof.pos]';
end
obj.elements = eleAll;

% intialize stiffness matrix and force vector
obj.K = zeros(nf_temp, nf_temp);
obj.Fp = zeros(nf_temp, 1);
obj.setPositions_F();

end