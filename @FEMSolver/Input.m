function Input(obj, runName)
%% read from input file
if exist('runName', 'var')
    run(runName);
else
    coordinates = [-1,-1; +1,-1; +1,+1; -1,+1];
    incidences = [1 1 2 3 4];
    eType = 'PhyElement_CPE4';
    EBC = [1, 1, 0;
        1, 2, 0;
        2, 2, 0;
        4, 1, 0];
    NBC = [2, 1, -15.0;
        3, 1, -15.0];
    nummat = 1;
    mults = 0.1*ones(10,1);
    mat = {'Hyperelastic',[40;60]};
end

%% check input file
nNodes_temp = size(coordinates, 1); % int, total nodes
ne_temp = size(incidences, 1); % int, total elements
ndof_temp = nNodes_temp * obj.dim; % int, total dofs
np_temp = size(EBC, 1); % int, total prescribed dofs
nf_temp = ndof_temp - np_temp; % int, total free dofs
EBC_temp = zeros(np_temp, 2);
NBC_temp = zeros(size(NBC, 1), 2);
FBC_temp = zeros(nf_temp, 1);

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
    % dofAll(currDof).v = EBC(temp, 3);
    dofAll(currDof).f = 0;
    dofAssigned(currDof) = 1;
    EBC_temp(temp, :) = [currDof, EBC(temp, 3)]; % collect essential BC
end

% initialize free dofs
temp2 = 1;
for temp = 1:ndof_temp
    if dofAssigned(temp) == 0
        dofAll(temp).p = 0;
        dofAll(temp).pos = temp2;
        dofAll(temp).v = 0;
        dofAll(temp).f = 0;
        FBC_temp(temp2) = temp; % collecting all natural BC
        temp2 = temp2 + 1;
    end
end
for temp = 1:size(NBC, 1)
    currDof = dofF(temp);
    dofAll(currDof).p = 0;
    dofAll(currDof).v = 0;
    % dofAll(currDof).f = NBC(temp, 3);
    NBC_temp(temp, :) = [currDof, NBC(temp, 3)];
end
% obj.dofs = dofAll(dofAssigned == 0);
obj.dofs = dofAll; % I fount that this way is easier for controlling stepped EBC and NBC
%
% rows have to be sorted so that assemble of residual term can be written
% in Matlab style
obj.NBC = sortrows(NBC_temp);
obj.EBC = sortrows(EBC_temp);
obj.FBC = sortrows(FBC_temp);
%
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
%
% initialize elements
eleAll(ne_temp, 1) = eval(eType); % avoid using 'switch'
for temp = 1:ne_temp
    eleAll(temp).id = incidences(temp, 1);
    eleAll(temp).eNodes = obj.nodes(incidences(temp, 2:end));
    eleAll(temp).nedof = reshape([eleAll(temp).eNodes.ndof], [], 1);
    eleAll(temp).dofMap = [eleAll(temp).nedof.pos]';
end
obj.elements = eleAll;
%
% intialize stiffness matrix and force vector
obj.K = zeros(nf_temp, nf_temp);
obj.Fp = zeros(nf_temp, 1);
obj.setPositions_F();
%
% multiplication factor for each step
if exist('mults', 'var')
    obj.mults = mults;
    obj.nstep = length(mults);
end
%
% read material model and parameters (allow only one material model for now)
% MateDecl = strcat(mat{1},'(',mat{2},')');
% obj.matAll = eval(MateDecl);
matString = sprintf('%.2f;' , mat{2});
matString = strcat(mat{1}, '([', matString(1:end-1), '])');
obj.matAll = eval(matString);
obj.matAll.Initialize(obj.dim);
nhardening = obj.matAll.nhardening;
%
% initialize the element subroutine
obj.elements.Initialize(nhardening);
end