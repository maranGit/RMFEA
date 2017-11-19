%============================midterm=======================================
% ========================axial loading====================================
% dimension = 2;
% coordinates = [0,0; +1,0; +1,+1; 0,+1];
% incidences = [1 1 2 3 4];
% eType = 'PhyElement_NLCPE4';
% EBC = [1, 1, 0;
%        1, 2, 0;
%        4, 1, 0;
%        4, 2, 0;
%        2, 1, 0.2;
%        3, 1, 0.2];
% NBC = [2, 2, 0;
%        3, 2, 0];
% nummat = 1;
% mults = 0.1*ones(10,1);
% mat = {'FiniteHyperElastic',[40;40]};
%
% ======================shear loading======================================
dimension = 2;
coordinates = [0,0; +1,0; +1,+1; 0,+1];
incidences = [1 1 2 3 4];
eType = 'PhyElement_NLCPE4';
EBC = [1, 1, 0;
       1, 2, 0;
       2, 1, 0;
       2, 2, 0;
       3, 1, 0.2;
       4, 1, 0.2];
NBC = [3, 2, 0;
       4, 2, 0];
nummat = 1;
mults = 0.1*ones(10,1);
mat = {'FiniteHyperElastic',[40;40]};