% This code was generated using ADiGator version 1.4
% ©2010-2014 Matthew J. Weinstein and Anil V. Rao
% ADiGator may be obtained at https://sourceforge.net/projects/adigator/ 
% Contact: mweinstein@ufl.edu
% Bugs/suggestions may be reported to the sourceforge forums
%                    DISCLAIMER
% ADiGator is a general-purpose software distributed under the GNU General
% Public License version 3.0. While the software is distributed with the
% hope that it will be useful, both the software and generated code are
% provided 'AS IS' with NO WARRANTIES OF ANY KIND and no merchantability
% or fitness for any purpose or application.

function C = hyperSensitiveFun_ADiGatorJac(z)
global ADiGator_hyperSensitiveFun_ADiGatorJac
if isempty(ADiGator_hyperSensitiveFun_ADiGatorJac); ADiGator_LoadData(); end
Gator1Data = ADiGator_hyperSensitiveFun_ADiGatorJac.hyperSensitiveFun_ADiGatorJac.Gator1Data;
% ADiGator Start Derivative Computations
global psStuff nstates ncontrols gravity 
%User Line: global
%User Line: % Extract the number of LGR points, the LGR differentiation matrix,
%User Line: % the LGR points, and the LGR weights
D = psStuff.D;
%User Line: D = psStuff.D;
tau = psStuff.tau;
%User Line: tau = psStuff.tau;
w = psStuff.w;
%User Line: w = psStuff.w;
%User Line: % Decompose the NLP decision vector into pieces containing
%User Line: %  - the state
%User Line: %  - the control
%User Line: %  - the initial time
%User Line: %  - the final time
cada1f1 = length(tau);
N.f = cada1f1 - 1;
%User Line: N = length(tau)-1;
cada1f1 = N.f + 1;
cada1f2 = nstates*cada1f1;
stateIndices.f = 1:cada1f2;
%User Line: stateIndices = 1:nstates*(N+1);
cada1f1 = N.f + 1;
cada1f2 = nstates*cada1f1;
cada1f3 = cada1f2 + 1;
cada1f4 = N.f + 1;
cada1f5 = nstates*cada1f4;
cada1f6 = ncontrols*N.f;
cada1f7 = cada1f5 + cada1f6;
controlIndices.f = cada1f3:cada1f7;
%User Line: controlIndices = (nstates*(N+1)+1):(nstates*(N+1)+ncontrols*N);
cada1f1 = length(controlIndices.f);
cada1f2 = controlIndices.f(cada1f1);
t0Index.f = cada1f2 + 1;
%User Line: t0Index = controlIndices(end)+1;
tfIndex.f = t0Index.f + 1;
%User Line: tfIndex = t0Index+1;
stateVector.dz0 = z.dz0(Gator1Data.Index1);
stateVector.f = z.f(stateIndices.f);
%User Line: stateVector = z(stateIndices);
controlVector.dz0 = z.dz0(Gator1Data.Index2);
controlVector.f = z.f(controlIndices.f);
%User Line: controlVector = z(controlIndices);
t0.dz0 = z.dz0(322);
t0.f = z.f(t0Index.f);
%User Line: t0 = z(t0Index);
tf.dz0 = z.dz0(323);
tf.f = z.f(tfIndex.f);
%User Line: tf = z(tfIndex);
%User Line: % Reshape the state and control parts of the NLP decision vector
%User Line: % to matrices of size (N+1) by nstates and N by ncontrols, respectively.
cada1f1 = N.f + 1;
state.dz0 = stateVector.dz0;
state.f = reshape(stateVector.f,cada1f1,nstates);
%User Line: state   = reshape(stateVector,N+1,nstates);
control.dz0 = controlVector.dz0;
control.f = reshape(controlVector.f,N.f,ncontrols);
%User Line: control = reshape(controlVector,N,ncontrols);
cada1f1 = size(state.f,1);
cada1f2 = cada1f1 - 1;
cada1f3 = 1:cada1f2;
stateLGR.dz0 = state.dz0(Gator1Data.Index3);
stateLGR.f = state.f(cada1f3,:);
%User Line: stateLGR = state(1:end-1,:);
%User Line: % Identify the different components of the state column-wise from stateLGR
x.dz0 = stateLGR.dz0(Gator1Data.Index4);
x.f = stateLGR.f(:,1);
%User Line: x = stateLGR(:,1);
u.dz0 = control.dz0; u.f = control.f;
%User Line: u = control;
%User Line: % Compute the right-hand side of the differential
%User Line: % equations at the LGR points
cada1f1dz0 = 3.*x.f(:).^(3-1).*x.dz0;
cada1f1 = x.f.^3;
cada1f2dz0 = -cada1f1dz0;
cada1f2 = uminus(cada1f1);
cada1td1 = zeros(320,1);
cada1td1(Gator1Data.Index5) = cada1f2dz0;
cada1td1(Gator1Data.Index6) = cada1td1(Gator1Data.Index6) + u.dz0;
rhs.dz0 = cada1td1;
rhs.f = cada1f2 + u.f;
%User Line: rhs = -x.^3+u;
%User Line: % Construct the defect constraints
cada1td1 = sparse(Gator1Data.Index7,Gator1Data.Index8,state.dz0,161,161);
cada1td1 = D*cada1td1;
cada1td1 = cada1td1(:);
cada1f1dz0 = full(cada1td1(Gator1Data.Index9));
cada1f1 = D*state.f;
cada1td1 = zeros(2,1);
cada1td1(2) = tf.dz0;
cada1td1(1) = cada1td1(1) + -t0.dz0;
cada1f2dz0 = cada1td1;
cada1f2 = tf.f - t0.f;
cada1tempdz0 = cada1f2dz0(Gator1Data.Index10);
cada1tf1 = rhs.f(Gator1Data.Index11);
cada1td1 = zeros(640,1);
cada1td1(Gator1Data.Index12) = cada1tf1(:).*cada1tempdz0;
cada1td1(Gator1Data.Index13) = cada1td1(Gator1Data.Index13) + cada1f2.*rhs.dz0;
cada1f3dz0 = cada1td1;
cada1f3 = cada1f2*rhs.f;
cada1f4dz0 = cada1f3dz0./2;
cada1f4 = cada1f3/2;
cada1td1 = zeros(1280,1);
cada1td1(Gator1Data.Index14) = cada1f1dz0;
cada1td1(Gator1Data.Index15) = cada1td1(Gator1Data.Index15) + -cada1f4dz0;
defects.dz0 = cada1td1;
defects.f = cada1f1 - cada1f4;
%User Line: defects = D*state - (tf-t0)*rhs/2;
%User Line: % Reshape the defect contraints into a column vector
cada1f1 = N.f*nstates;
defects.dz0 = defects.dz0;
defects.f = reshape(defects.f,cada1f1,1);
%User Line: defects = reshape(defects,N*nstates,1);
cada1td1 = zeros(2,1);
cada1td1(2) = tf.dz0;
cada1td1(1) = cada1td1(1) + -t0.dz0;
cada1f1dz0 = cada1td1;
cada1f1 = tf.f - t0.f;
cada1f2dz0 = cada1f1dz0./2;
cada1f2 = cada1f1/2;
cada1f3dz0 = 0.5.*cada1f2dz0;
cada1f3 = cada1f2*0.5;
cada1td1 = x.f(:).*x.dz0;
cada1td1 = cada1td1 + x.f(:).*x.dz0;
cada1f4dz0 = cada1td1;
cada1f4 = x.f.*x.f;
cada1td1 = u.f(:).*u.dz0;
cada1td1 = cada1td1 + u.f(:).*u.dz0;
cada1f5dz0 = cada1td1;
cada1f5 = u.f.*u.f;
cada1td1 = zeros(320,1);
cada1td1(Gator1Data.Index16) = cada1f4dz0;
cada1td1(Gator1Data.Index17) = cada1td1(Gator1Data.Index17) + cada1f5dz0;
cada1f6dz0 = cada1td1;
cada1f6 = cada1f4 + cada1f5;
cada1tempdz0 = cada1f3dz0(Gator1Data.Index18);
cada1tf1 = cada1f6(Gator1Data.Index19);
cada1td1 = zeros(640,1);
cada1td1(Gator1Data.Index20) = cada1tf1(:).*cada1tempdz0;
cada1td1(Gator1Data.Index21) = cada1td1(Gator1Data.Index21) + cada1f3.*cada1f6dz0;
Lagrangian.dz0 = cada1td1;
Lagrangian.f = cada1f3*cada1f6;
%User Line: Lagrangian = ((tf-t0)/2)*0.5*(x.*x+u.*u);
cada1f1 = psStuff.w.';
cada1td1 = sparse(Gator1Data.Index22,Gator1Data.Index23,Lagrangian.dz0,160,322);
cada1td1 = cada1f1*cada1td1;
cada1td1 = cada1td1(:);
J.dz0 = full(cada1td1(Gator1Data.Index24));
J.f = cada1f1*Lagrangian.f;
%User Line: J = psStuff.w.'*Lagrangian;
%User Line: % Construct the vector of constraints plus the objective function
cada1td1 = zeros(1602,1);
cada1td1(Gator1Data.Index25) = J.dz0;
cada1td1(Gator1Data.Index26) = defects.dz0;
C.dz0 = cada1td1;
C.f = [J.f;defects.f];
%User Line: C = [J; defects];
C.dz0_size = [161,323];
C.dz0_location = Gator1Data.Index27;
end


function ADiGator_LoadData()
global ADiGator_hyperSensitiveFun_ADiGatorJac
ADiGator_hyperSensitiveFun_ADiGatorJac = load('hyperSensitiveFun_ADiGatorJac.mat');
return
end