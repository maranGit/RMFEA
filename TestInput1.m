%=====================square element=======================================
% dimension = 2;
% coordinates = [-1,-1; +1,-1; +1,+1; -1,+1];
% incidences = [1 1 2 3 4];
% eType = 'PhyElement_CPE4';
% EBC = [1, 1, 0;
%        1, 2, 0;
%        4, 1, 0
%        4, 2, 0];
% NBC = [2, 2, 1;
%        3, 2, 1];
%=====================rectangular elements=================================
% dimension = 2;
% coordinates = [0,0; 1.2,0; 0,1; 1,1; 2,0; 2,1];
% incidences = [1 1 2 4 3; 2 2 5 6 4];
% eType = 'PhyElement_CPE4';
% EBC = [1, 1, 0;
%        1, 2, 0;
%        3, 1, 0;
%        2, 2, 0;
%        5, 2, 0];
% NBC = [5, 1, 1;
%        6, 1, 1];
%==============compare element stiffness with analytical result============
% dimension = 2;
% coordinates = [1,0; 2,0; 2.25,1.5; 1.25,1];
% incidences = [1 1 2 3 4];
% eType = 'PhyElement_CPE4';
% EBC = [1, 1, 0;
%        1, 2, 0;
%        3, 1, 0];
% NBC = [3, 2, 1];
%============================pure shear====================================
% dimension = 2;
% coordinates = [0,0; 1,0; 0,1; 1,1];
% incidences = [1 1 2 4 3];
% eType = 'PhyElement_CPE4';
% EBC = [1, 1, 0;
%        1, 2, 0;
%        2, 2, 0];
% NBC = [2, 1, -1;
%        3, 1, 1;
%        3, 2, -1;
%        4, 1, 1;
%        4, 2, 1];
%============================patch test====================================
% dimension = 2;
% coordinates = [0, 0;
%                0.7, 0;
%                2, 0;
%                0, 1;
%                1.2, 0.7;
%                2, 1;
%                0, 2;
%                1.5, 2;
%                2, 2];
% incidences = [1 1 2 5 4;
%               2 2 3 6 5;
%               3 4 5 8 7;
%               4 5 6 9 8];
% eType = 'PhyElement_CPE4';
% EBC = [1, 1, 0;
%        1, 2, 0;
%        4, 1, 0];
% NBC = [3, 1, 1;
%        6, 1, 2;
%        7, 1, -1;
%        9, 1, 1];
%============================HW1 Q5========================================
dimension = 2;
coordinates = [-1,-1; +1,-1; +1,+1; -1,+1];
incidences = [1 1 2 3 4];
eType = 'PhyElement_CPE4';
EBC = [1, 1, 0;
       1, 2, 0;
       4, 1, 0;
       2, 2, 0;
       2, 1, 0.02;
       3, 1, 0.02];
NBC = [4, 2, 0;
       3, 2, 0];