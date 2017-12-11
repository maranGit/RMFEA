strainTOL =100;1000;1;
% =========================================================================
if strainTOL == 100
% single element test
dimension = 2;
coordinates = [-1,-1; +1,-1; +1,+1; -1,+1];
incidences = [1 1 2 3 4];
eType = 'Hypo2d';
EBC = [1, 1, 0;
       1, 2, 0;
       2, 2, 0;
       4, 1, 0;
       2, 1, 0.4;
       3, 1, 0.4];
NBC = [3, 2, 0;
       4, 2, 0];
lint = 4;
nummat = 1;
mults = 0.02*ones(50,1);
mat = {'VonMises_finite',[12000;0.25;22;2.5;3.5]};
% mat = {'Hypoelastic',[12000;0.3]};
% =========================================================================
elseif strainTOL == 1
% small strain compare
dimension = 2;
coordinates = [-1,-1; +1,-1; +1,+1; -1,+1];
incidences = [1 1 2 3 4];
eType = 'PhyElement_CPE4';
EBC = [1, 1, 0;
       1, 2, 0;
       2, 2, 0;
       4, 1, 0;
       2, 1, 0.04;
       3, 1, 0.04];
NBC = [3, 2, 0;
       4, 2, 0];
lint = 4;
nummat = 1;
mults = 0.1*ones(10,1);
% mat = {'VonMises',[1200;0.25;22;2.5;3.5]};
mat = {'VonMises',[12000;0.3;10000;1000;3]};
% =========================================================================
elseif strainTOL == 1000
% multi elements test
dimension = 2;
coordinates = [0, 0; 1, 0; 2, 0; 3, 0;
    0, 1; 1, 1; 2, 1; 3, 1;
    0, 2; 1, 2; 2, 2; 3, 2;
    0, 3; 1, 3; 2, 3; 3, 3];
incidences = [1 0 1 5 4;
              2 1 2 6 5;
              3 2 3 7 6;
              4 4 5 9 8;
              5 5 6 10 9;
              6 6 7 11 10;
              7 8 9 13 12;
              8 9 10 14 13;
              9 10 11 15 14];
incidences(:, 2:5) = incidences(:, 2:5) + 1;
eType = 'Hypo2d';
% eType = 'ele3';
EBC = [1, 1, 0;
       1, 2, 0;
       5, 1, 0;
       9, 1, 0;
       13, 1, 0;
       2, 2, 0;
       3, 2, 0;
       4, 2, 0;
       4, 1, 0.6;
       8, 1, 0.6;
       12, 1, 0.6;
       16, 1, 0.6];
NBC = [13, 2, 0;
       14, 2, 0;
       15, 2, 0;
       16, 2, 0];
lint = 4;
nummat = 1;
mults = 0.02*ones(50,1);
mat = {'VonMises_finite',[12000;0.3;100;1000;3]};
% mat = {'Hypoelastic',[1200;0.3]};
end