function d = Solve_Dofs(obj)
d = obj.K \ obj.F;
end