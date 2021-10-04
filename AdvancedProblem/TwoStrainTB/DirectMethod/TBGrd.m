function grd = TBGrd(Z)
% computes the gradient

output = TBObj_Jac(Z);
grd    = output;

end

