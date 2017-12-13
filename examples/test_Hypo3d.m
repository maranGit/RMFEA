% single element test
dimension = 3;
coordinates = [0,0,0; 1,0,0; 1,1,0; 0,1,0; 0,0,1; 1,0,1; 1,1,1; 0,1,1];
incidences = [1 1 2 3 4 5 6 7 8];
eType = 'Hypo3d';
EBC = [1, 1, 0;
       1, 2, 0;
       1, 3, 0;
       4, 1, 0;
       5, 1, 0;
       8, 1, 0;
       2, 2, 0;
       5, 2, 0;
       6, 2, 0;
       2, 3, 0;
       3, 3, 0;
       4, 3, 0;
       5, 3, 0;
       6, 3, 0;
       7, 3, 0;
       8, 3, 0;
       2, 1, 0.2;
       3, 1, 0.2;
       6, 1, 0.2;
       7, 1, 0.2];
NBC = [3, 2, 0;
       4, 2, 0;
       7, 2, 0;
       8, 2, 0];
lint = 8;
nummat = 1;
mults = 0.02*ones(50,1);
mat = {'VonMises_finite',[12000;0.25;22;2.5;3.5]};