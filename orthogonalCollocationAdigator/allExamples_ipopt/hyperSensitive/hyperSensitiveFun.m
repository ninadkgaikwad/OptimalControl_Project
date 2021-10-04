function C = hyperSensitiveFun(z)

global psStuff nstates ncontrols 

% Extract the number of LGR points, the LGR differentiation matrix, 
% the LGR points, and the LGR weights
D = psStuff.D;
tau = psStuff.tau;
w = psStuff.w;

% Decompose the NLP decision vector into pieces containing
%  - the state
%  - the control
%  - the initial time
%  - the final time
N = length(tau)-1;
stateIndices = 1:nstates*(N+1);
controlIndices = (nstates*(N+1)+1):(nstates*(N+1)+ncontrols*N);
t0Index = controlIndices(end)+1;
tfIndex = t0Index+1;
stateVector = z(stateIndices);
controlVector = z(controlIndices);
t0 = z(t0Index);
tf = z(tfIndex);

% Reshape the state and control parts of the NLP decision vector
% to matrices of size (N+1) by nstates and N by ncontrols, respectively.
state   = reshape(stateVector,N+1,nstates);
control = reshape(controlVector,N,ncontrols);
stateLGR = state(1:end-1,:);

% Identify the different components of the state column-wise from stateLGR
x = stateLGR(:,1);
u = control;

% Compute the right-hand side of the differential 
% equations at the LGR points
rhs = -x.^3+u;

% Construct the defect constraints
defects = D*state - (tf-t0)*rhs/2;

% Reshape the defect contraints into a column vector
defects = reshape(defects,N*nstates,1);
Lagrangian = ((tf-t0)/2)*0.5*(x.*x+u.*u);
J = psStuff.w.'*Lagrangian;

% Construct the vector of constraints plus the objective function
C = [J; defects];
