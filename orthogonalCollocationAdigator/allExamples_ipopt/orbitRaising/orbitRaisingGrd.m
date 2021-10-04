function grd = orbitRaisingGrd(Z)
% computes the gradient

output = orbitRaisingObj_Jac(Z);
grd    = output;

end

