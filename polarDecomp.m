function [R, U] = polarDecomp(F)
if ~isequal(size(F), [3,3])
    error('F should be a 3x3 matrix');
end
if det(F) < 1e-6
    error('determinant of F should be positive');
end
C = transpose(F)*F;
lamda2 = eig(C);
lamda = sqrt(lamda2);
i1 = sum(lamda);
i2 = lamda(1)*lamda(2)+lamda(1)*lamda(3)+lamda(3)*lamda(2);
i3 = lamda(1)*lamda(2)*lamda(3);
D = i1*i2-i3;
U = (-C*C+(i1*i1-i2)*C+i1*i3*eye(3))/D;
Uinv = (C-i1*U+i2*eye(3))/i3;
R=F*Uinv;
end