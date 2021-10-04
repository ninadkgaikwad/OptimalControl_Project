function [C,J] = orbitRaisingFunWrapper(z);

% zad = gradientinit(z);
% F = orbitRaisingFun(zad);
% C = F.x;
% J = F.dx;

global iGfun jGvar

[J,C] = orbitRaisingFun_Jac(z);
G     = snfindG(iGfun,jGvar,J);