function C = TBFun(z)

%-----------------------------------------------------------------%
% Objective and constraint functions for the orbit-raising        %
% problem.  This function is designed to be used with the NLP     %
% solver SNOPT.                                                   %
%-----------------------------------------------------------------%
%      DO NOT FOR ANY REASON ALTER THE LINE OF CODE BELOW!        %
global psStuff nstates ncontrols npaths 
global Beta1 Beta2 Beta3 Mu
global d1 d2 k1 k2 r1 r2 
global p q N1 B1 B2 
%      DO NOT FOR ANY REASON ALTER THE LINE OF CODE ABOVE!        %
%-----------------------------------------------------------------%

%-----------------------------------------------------------------%
% Radau pseudospectral method quantities required:                %
%   - Differentiation matrix (psStuff.D)                          %
%   - Legendre-Gauss-Radau weights (psStuff.w)                    %
%   - Legendre-Gauss-Radau points (psStuff.tau)                   %
%-----------------------------------------------------------------%
D = psStuff.D; tau = psStuff.tau; w = psStuff.w;

%-----------------------------------------------------------------%
% Decompose the NLP decision vector into pieces containing        %
%    - the state                                                  %
%    - the control                                                %
%    - the initial time                                           %
%    - the final time                                             %
%-----------------------------------------------------------------%
N = length(tau)-1;
stateIndices = 1:nstates*(N+1);
controlIndices = (nstates*(N+1)+1):(nstates*(N+1)+ncontrols*N);
t0Index = controlIndices(end)+1;
tfIndex = t0Index+1;
stateVector = z(stateIndices);
controlVector = z(controlIndices);
t0 = z(t0Index);
tf = z(tfIndex);

%-----------------------------------------------------------------%
% Reshape the state and control parts of the NLP decision vector  %
% to matrices of sizes (N+1) by nstates and (N+1) by ncontrols,   %
% respectively.  The state is approximated at the N LGR points    %
% plus the final point.  Thus, each column of the state vector is %
% length N+1.  The LEFT-HAND SIDE of the defect constraints, D*X, %
% uses the state at all of the points (N LGR points plus final    %
% point).  The RIGHT-HAND SIDE of the defect constraints,         %
% (tf-t0)F/2, uses the state and control at only the LGR points.  %
% Thus, it is necessary to extract the state approximations at    %
% only the N LGR points.  Finally, in the Radau pseudospectral    %
% method, the control is approximated at only the N LGR points.   %
%-----------------------------------------------------------------%
statePlusEnd   = reshape(stateVector,N+1,nstates);
control = reshape(controlVector,N,ncontrols);
stateLGR = statePlusEnd(1:end-1,:);

%-----------------------------------------------------------------%
% Identify the components of the state column-wise from stateLGR. % 
%-----------------------------------------------------------------%
S = stateLGR(:,1);
T = stateLGR(:,2);
L1 = stateLGR(:,3);
L2 = stateLGR(:,4);
I1 = stateLGR(:,5);
I2 = stateLGR(:,6);
u1 = control(:,1);
u2 = control(:,2);

%-----------------------------------------------------------------%
% Compute the right-hand side of the differential equations at    %
% the N LGR points.  Each component of the right-hand side is     %
% stored as a column vector of length N, that is each column has  %
% the form                                                        %
%                   [ f_i(x_1,u_1,t_1) ]                          %
%                   [ f_i(x_2,u_2,t_2) ]                          %
%                           .                                     %
%                           .                                     %
%                           .                                     %
%                   [ f_i(x_N,u_N,t_N) ]                          %
% where "i" is the right-hand side of the ith component of the    %
% vector field f.  It is noted that in MATLABB the calculation of %
% the right-hand side is vectorized.                              %
%-----------------------------------------------------------------%

% Computing Derivatives
S_Derivative = (Mu*N1) - (Beta1*S.*(I1/N1)) - (Beta3*S.*(I2/N1)) - (Mu*S);

T_Derivative = (r1*u1.*L1) - (Mu*T) - (r2*(p+q)*(1-(1-u2)).*I1)...
                - (Beta2*T.*(I1/N1)) - (Beta3*T.*(I2/N1));

L1_Derivative = (Beta1*S.*(I1/N1)) - ((Mu+k1)*L1) - (r1*u1.*L1)...
                 + (p*r2*(1-u2).*I1) + (Beta2*T.*(I1/N1))...
                 - (Beta3*L1.*(I2/N1));

L2_Derivative = (q*r2*(1-u2).*I1) - ((Mu+k2)*L2) ...
                + (Beta3*(S+L1+T).*(I2/N1));

I1_Derivative = (k1*L1) - ((Mu+d1)*I1) - (r2*I1);

I2_Derivative = (k2*L2) - ((Mu+d2)*I2);


diffeqRHS = [S_Derivative, T_Derivative, L1_Derivative, ...
    L2_Derivative, I1_Derivative, I2_Derivative];

%-----------------------------------------------------------------%
% Compute the left-hand side of the defect constraints, recalling %
% that the left-hand side is computed using the state at the LGR  %
% points PLUS the final point.                                    %
%-----------------------------------------------------------------%
diffeqLHS = D*statePlusEnd;

%-----------------------------------------------------------------%
% Construct the defect constraints at the N LGR points.           %
% Remember that the right-hand side needs to be scaled by the     %
% factor (tf-t0)/2 because the rate of change of the state is     %
% being taken with respect to $\tau\in[-1,+1]$.  Thus, we have    %
% $dt/t\dau=(tf-t0)/2$.                                           %
%-----------------------------------------------------------------%
defects = diffeqLHS-(tf-t0)*diffeqRHS/2;

%-----------------------------------------------------------------%
% Construct the path constraints at the N LGR points.             %
%-----------------------------------------------------------------%
paths = S + T + L1 + L2 + I1 + I2 - N1;

%-----------------------------------------------------------------%
% Reshape the defect contraints into a column vector.             % 
%-----------------------------------------------------------------%
defects = reshape(defects,N*nstates,1);

%-----------------------------------------------------------------%
% Reshape the Path contraints into a column vector.             % 
%-----------------------------------------------------------------%
paths = reshape(paths,N*npaths,1);

%-----------------------------------------------------------------%
% Construct the objective function plus constraint vector.        %
%-----------------------------------------------------------------%
J   = sum(L2 + I2 + ((1/2)*B1*u1.^2) + ((1/2)*B2*u2.^2));

C = [J; defects; paths];
