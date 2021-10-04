function obj = orbitRaisingObj(z)
% Computes the objective function of the problem

global psStuff nstates ncontrols CONSTANTS       

%-----------------------------------------------------------------%
%         Extract the constants used in the problem.              %
%-----------------------------------------------------------------%
MU = CONSTANTS.MU; mdot = CONSTANTS.mdot; T = CONSTANTS.T;

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
t = (tf-t0)*(tau+1)/2+t0;
tLGR = t(1:end-1);

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
stateLGR = statePlusEnd(1:end-1,:);
control = reshape(controlVector,N,ncontrols);

%-----------------------------------------------------------------%
% Identify the components of the state column-wise from stateLGR. % 
%-----------------------------------------------------------------%
r = stateLGR(:,1);
theta = stateLGR(:,1);
vr = stateLGR(:,3);
vtheta = stateLGR(:,4);
m = stateLGR(:,5);
u1 = control(:,1);
u2 = control(:,2);

%-----------------------------------------------------------------%
% The quantity STATEF is the value of the state at the final      %
% time, tf, which corresponds to the state at $\tau=1$.           %
%-----------------------------------------------------------------%
stateF = statePlusEnd(end,:);
%-----------------------------------------------------------------%
% The orbit-raising problem contains one nonlinear boundary       %
% condition $\sqrt{mu/r(t_f)-v_\theta(t_f) = 0$.  Because $r(t)$  %
% and $v_\theta(t)$ are the first and fourth components of the    %
% state, it is necessary to extract stateF(1) and stateF(4) in    %
% order to compute this boundary condition function.              %
%-----------------------------------------------------------------%
rF = stateF(1);

% Cost Function
J   = -rF;
obj = J;

end

