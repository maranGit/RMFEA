% Ran Ma
% FEA

clear all
clc
fclose('all');
%
% copied from NLFEA
if exist('filename.mat','file')
    clear load
    load('filename.mat', 'defaultname','defaultpath')
else
    defaultname = 'DNE';
end
fprintf('Provide name of input file (default is %s)\n',defaultname);
[filename,pathname] = uigetfile('*.m','Select the NLFEA input file');
if isnumeric(filename)
    pathname = defaultpath;
    filename = defaultname;
end
defaultname = filename;
defaultpath = pathname;
save('filename.mat','defaultname','defaultpath');
run([pathname filename(1:end-2)])
load('filename.mat') % reload filename so that the workspace knows which file it was
% copied from NLFEA
%
runName = filename;
if (~exist('runName','var'))
    runName = 'TestInput1';
end
femSolver = FEMSolver();
femSolver.FEMSolve(runName);