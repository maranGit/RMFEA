% single element test
dimension = 3;
coordinates = [0,0,0; 1,0,0; 1,1,0; 0,1,0; 0,0,1; 1,0,1; 1,1,1; 0,1,1];
incidences = [1 1 2 3 4 5 6 7 8];
eType = 'Small3d';
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
       2, 1, 0.02;
       3, 1, 0.02;
       6, 1, 0.02;
       7, 1, 0.02;];
NBC = [3, 2, 0;
       4, 2, 0;
       7, 2, 0;
       8, 2, 0];
lint = 8;
nummat = 1;
mults = 0.01*ones(100,1);
mat = {'VonMises',[12000;0.3;100;1000;3]};
%%
% multi elements test
% dimension = 3;
% coordinates = [0, 0; 1, 0; 2, 0; 3, 0;
%     0, 1; 1, 1; 2, 1; 3, 1;
%     0, 2; 1, 2; 2, 2; 3, 2;
%     0, 3; 1, 3; 2, 3; 3, 3];
% incidences = [1 0 1 5 4;
%               2 1 2 6 5;
%               3 2 3 7 6;
%               4 4 5 9 8;
%               5 5 6 10 9;
%               6 6 7 11 10;
%               7 8 9 13 12;
%               8 9 10 14 13;
%               9 10 11 15 14];
% incidences(:, 2:5) = incidences(:, 2:5) + 1;
% eType = 'Small2d';
% EBC = [1, 1, 0;
%        1, 2, 0;
%        5, 1, 0;
%        9, 1, 0;
%        13, 1, 0;
%        2, 2, 0;
%        3, 2, 0;
%        4, 2, 0;
%        4, 1, 0.2;
%        8, 1, 0.2;
%        12, 1, 0.2;
%        16, 1, 0.2];
% NBC = [13, 2, 0;
%        14, 2, 0;
%        15, 2, 0;
%        16, 2, 0];
% lint = 4;
% nummat = 1;
% mults = [0.02*ones(50,1)];
% mat = {'VonMises',[12000;0.3;100;1000;3]};