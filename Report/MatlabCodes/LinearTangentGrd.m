function grd = LinearTangentGrd(Z)
% computes the gradient

output = LinearTangentObj_Jac(Z);
grd    = output;

end

