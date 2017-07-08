dimension = 2;
coordinates = [0,0; 1,0; 0,1; 1,1; 2,0; 2,1];
incidences = [1 1 2 4 3; 2 2 5 6 4];
eType = 'PhyElement_CPE4';
EBC = [1, 1, 0;
       1, 2, 0;
       3, 1, 0];
NBC = [5, 1, 1;
       6, 1, 1];
