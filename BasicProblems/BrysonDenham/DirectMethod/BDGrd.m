function grd = BDGrd(Z)
% computes the gradient

output = BDObj_Jac(Z);
grd    = output;

end

