function grd = RobotArmGrd(Z)
% computes the gradient

output = RobotArmObj_Jac(Z);
grd    = output;

end

