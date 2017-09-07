# RMFEA
Illustration for Homework 1, Problem 5
1. Essential BC and natural BC are in "./TestInput1.m";
2. Material parameters are hardcoded in "./@PhyElement_CPE4/PhyElement_CPE4.m" now. They will be moved to material module when more complicated model is implemented;
3. To run the simulation: 
	test = FEMSolver(2);
	test.FEMSolve('TestInput1.m');
4. check stress at integration point:
	test.elements.stress;
	format:
	xx1, xx2, xx3, xx4
	yy1, yy2, yy3, yy4
	xy1, xy2, xy3, xy4