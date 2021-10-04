function grd = hyperSensitiveGrd(Z)
% computes the gradient

output = hyperSensitiveObj_Jac(Z);
grd    = output;

end

