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

function obj = RobotArmObj_ADiGatorJac(z)
global ADiGator_RobotArmObj_ADiGatorJac
if isempty(ADiGator_RobotArmObj_ADiGatorJac); ADiGator_LoadData(); end
Gator1Data = ADiGator_RobotArmObj_ADiGatorJac.RobotArmObj_ADiGatorJac.Gator1Data;
% ADiGator Start Derivative Computations
%User Line: % Computes the objective function of the problem
global psStuff nstates ncontrols L 
%User Line: global
%User Line: %-----------------------------------------------------------------%
%User Line: % Radau pseudospectral method quantities required:                %
%User Line: %   - Differentiation matrix (psStuff.D)                          %
%User Line: %   - Legendre-Gauss-Radau weights (psStuff.w)                    %
%User Line: %   - Legendre-Gauss-Radau points (psStuff.tau)                   %
%User Line: %-----------------------------------------------------------------%
D = psStuff.D;
%User Line: D = psStuff.D;
tau = psStuff.tau;
%User Line: tau = psStuff.tau;
w = psStuff.w;
%User Line: w = psStuff.w;
%User Line: %-----------------------------------------------------------------%
%User Line: % Decompose the NLP decision vector into pieces containing        %
%User Line: %    - the state                                                  %
%User Line: %    - the control                                                %
%User Line: %    - the initial time                                           %
%User Line: %    - the final time                                             %
%User Line: %-----------------------------------------------------------------%
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
t0.dz0 = z.dz0(367);
t0.f = z.f(t0Index.f);
%User Line: t0 = z(t0Index);
tf.dz0 = z.dz0(368);
tf.f = z.f(tfIndex.f);
%User Line: tf = z(tfIndex);
%User Line: %-----------------------------------------------------------------%
%User Line: % Reshape the state and control parts of the NLP decision vector  %
%User Line: % to matrices of sizes (N+1) by nstates and (N+1) by ncontrols,   %
%User Line: % respectively.  The state is approximated at the N LGR points    %
%User Line: % plus the final point.  Thus, each column of the state vector is %
%User Line: % length N+1.  The LEFT-HAND SIDE of the defect constraints, D*X, %
%User Line: % uses the state at all of the points (N LGR points plus final    %
%User Line: % point).  The RIGHT-HAND SIDE of the defect constraints,         %
%User Line: % (tf-t0)F/2, uses the state and control at only the LGR points.  %
%User Line: % Thus, it is necessary to extract the state approximations at    %
%User Line: % only the N LGR points.  Finally, in the Radau pseudospectral    %
%User Line: % method, the control is approximated at only the N LGR points.   %
%User Line: %-----------------------------------------------------------------%
cada1f1 = N.f + 1;
statePlusEnd.dz0 = stateVector.dz0;
statePlusEnd.f = reshape(stateVector.f,cada1f1,nstates);
%User Line: statePlusEnd   = reshape(stateVector,N+1,nstates);
control.dz0 = controlVector.dz0;
control.f = reshape(controlVector.f,N.f,ncontrols);
%User Line: control = reshape(controlVector,N,ncontrols);
cada1f1 = size(statePlusEnd.f,1);
cada1f2 = cada1f1 - 1;
cada1f3 = 1:cada1f2;
stateLGR.dz0 = statePlusEnd.dz0(Gator1Data.Index3);
stateLGR.f = statePlusEnd.f(cada1f3,:);
%User Line: stateLGR = statePlusEnd(1:end-1,:);
%User Line: %-----------------------------------------------------------------%
%User Line: % Identify the components of the state column-wise from stateLGR. %
%User Line: %-----------------------------------------------------------------%
x1.dz0 = stateLGR.dz0(Gator1Data.Index4);
x1.f = stateLGR.f(:,1);
%User Line: x1 = stateLGR(:,1);
x2.dz0 = stateLGR.dz0(Gator1Data.Index5);
x2.f = stateLGR.f(:,2);
%User Line: x2 = stateLGR(:,2);
x3.dz0 = stateLGR.dz0(Gator1Data.Index6);
x3.f = stateLGR.f(:,3);
%User Line: x3 = stateLGR(:,3);
x4.dz0 = stateLGR.dz0(Gator1Data.Index7);
x4.f = stateLGR.f(:,4);
%User Line: x4 = stateLGR(:,4);
x5.dz0 = stateLGR.dz0(Gator1Data.Index8);
x5.f = stateLGR.f(:,5);
%User Line: x5 = stateLGR(:,5);
x6.dz0 = stateLGR.dz0(Gator1Data.Index9);
x6.f = stateLGR.f(:,6);
%User Line: x6 = stateLGR(:,6);
u1.dz0 = control.dz0(Gator1Data.Index10);
u1.f = control.f(:,1);
%User Line: u1 = control(:,1);
u2.dz0 = control.dz0(Gator1Data.Index11);
u2.f = control.f(:,2);
%User Line: u2 = control(:,2);
u3.dz0 = control.dz0(Gator1Data.Index12);
u3.f = control.f(:,3);
%User Line: u3 = control(:,3);
%User Line: % Cost function
J.dz0 = tf.dz0; J.f = tf.f;
%User Line: J   = tf;
obj.dz0 = J.dz0; obj.f = J.f;
%User Line: obj = J;
obj.dz0_size = 368;
obj.dz0_location = Gator1Data.Index13;
end


function ADiGator_LoadData()
global ADiGator_RobotArmObj_ADiGatorJac
ADiGator_RobotArmObj_ADiGatorJac = load('RobotArmObj_ADiGatorJac.mat');
return
end