function [C,G] = brachistochroneFunWrapper(z);

global iGfun jGvar

[J,C] = brachistochroneFun_Jac(z);
G     = snfindG(iGfun,jGvar,J);
