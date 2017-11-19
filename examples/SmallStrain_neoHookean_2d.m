dimension = 2;
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