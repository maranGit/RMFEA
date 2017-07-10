function Calculate_Stress_Strain(obj)
E1 = obj.E1;
E2 = obj.E2;
G = obj.G;

% material stiffness
D = [E1, E2, 0;
     E2, E1, 0;
     0,  0,  G];
 
disp = transpose( [obj.nedof.v] );
obj.fde = obj.ke * disp;
obj.strain = [obj.Bmat_1 * disp, obj.Bmat_2 * disp, ...
    obj.Bmat_3 * disp, obj.Bmat_4 * disp];
obj.stress = D * obj.strain;
end