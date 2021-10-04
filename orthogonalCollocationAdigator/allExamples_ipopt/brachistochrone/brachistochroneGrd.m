function grd = brachistochroneGrd(Z)
% computes the gradient

output = brachistochroneObj_Jac(Z);
grd    = output;

end

