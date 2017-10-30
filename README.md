# RMFEA
Illustration for Midterm:
1. type 'RMFEA' in command window
2. choose 'MidTermInput.m' as input file. Change MidTermInput.m to transfer between axial loading and shear loading
3. check output file in 'Output' for results collection
4. major revised files: (1) 'Calculate_ElementStiffness_Force.m' in folder '@PhyElement_NLCPE4'
                        (2) 'FiniteHyperElastic.m' in folder 'RMFEA'
5. results verification:
   (1) displacement = [femSolver.dofs.v];
   (2) CauchyStress = [femSolver.elements.stress]; % stress at gp1, gp2, gp3, gp4 in Voigt notation.