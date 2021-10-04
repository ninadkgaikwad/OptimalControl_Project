function [C,G] = hyperSensitiveFunWrapper(z);

global iGfun jGvar

[J,C] = hyperSensitiveFun_Jac(z);
G     = snfindG(iGfun,jGvar,J);
