function obj = TBObj(z)
% Computes the objective function of the problem

global psStuff nstates ncontrols 
global Beta1 Beta2 Beta3 Mu
global d1 d2 k1 k2 r1 r2 
global p q N1 B1 B2 

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

% Cost function
J   = sum(L2 + I2 + ((1/2)*B1*u1.^2) + ((1/2)*B2*u2.^2));
obj = J;

end

