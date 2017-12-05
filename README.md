# RMFEA
Illustration for HW5:
1. type 'RMFEA' in command window
2. choose 'smallStrain_J2.m' in 'examples' folder as input file.
3. check output file in 'Output' for results collection
4. major revised files: (1) 'Calculate_ElementStiffness_Force.m' in folder '@PhyElement_CPE4', add internal variables
                        (2) 'VonMises.m' in folder 'material'
5. results verification:
   (1) displacement = [femSolver.dofs.v];
   (2) CauchyStress = [femSolver.elements.stress]; % stress at gp1, gp2, gp3, gp4 in Voigt notation.