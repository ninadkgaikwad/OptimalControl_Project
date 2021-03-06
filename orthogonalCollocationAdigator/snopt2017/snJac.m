function  [A,iAfun,jAvar,iGfun,jGvar] = snJac(userfun,x0,xlow,xupp,nF)
%function [A,iAfun,jAvar,iGfun,jGvar] = snJac(usrfun,x0,xlow,xupp,nF)
%         Finds the coordinate structure for the Jacobian.

findJacOption = 17;

userFG = checkFun(userfun,'SNOPT');

if nargout(userFG) ~= 1,
  error('SNOPT:InputArgs','%s should only compute the objective function',inputname(userfun));
end

[A,iAfun,jAvar,iGfun,jGvar] = snoptmex(findJacOption,userFG,x0,xlow,xupp,nF);
