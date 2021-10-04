function grd = GMRGrd(Z)
% computes the gradient

output = GMRObj_Jac(Z);
grd    = output;

end

