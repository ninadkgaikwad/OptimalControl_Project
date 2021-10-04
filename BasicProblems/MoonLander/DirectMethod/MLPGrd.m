function grd = MLPGrd(Z)
% computes the gradient

output = MLPObj_Jac(Z);
grd    = output;

end

