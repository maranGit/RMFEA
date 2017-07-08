% Ran Ma
% FEA

if (~exist('runName','var'))
    runName = 'TestInput1';
end
femSolver = FEMSolver(2);
femSolver.FEMSolve(runName);