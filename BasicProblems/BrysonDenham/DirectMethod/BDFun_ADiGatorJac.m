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

function C = BDFun_ADiGatorJac(z)
global ADiGator_BDFun_ADiGatorJac
if isempty(ADiGator_BDFun_ADiGatorJac); ADiGator_LoadData(); end
Gator1Data = ADiGator_BDFun_ADiGatorJac.BDFun_ADiGatorJac.Gator1Data;
% ADiGator Start Derivative Computations
%User Line: %-----------------------------------------------------------------%
%User Line: % Objective and constraint functions for the orbit-raising        %
%User Line: % problem.  This function is designed to be used with the NLP     %
%User Line: % solver SNOPT.                                                   %
%User Line: %-----------------------------------------------------------------%
%User Line: %      DO NOT FOR ANY REASON ALTER THE LINE OF CODE BELOW!        %
global psStuff nstates ncontrols l 
%User Line: global
%User Line: %      DO NOT FOR ANY REASON ALTER THE LINE OF CODE ABOVE!        %
%User Line: %-----------------------------------------------------------------%
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
t0.dz0 = z.dz0(123);
t0.f = z.f(t0Index.f);
%User Line: t0 = z(t0Index);
tf.dz0 = z.dz0(124);
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
x.dz0 = stateLGR.dz0(Gator1Data.Index4);
x.f = stateLGR.f(:,1);
%User Line: x = stateLGR(:,1);
v.dz0 = stateLGR.dz0(Gator1Data.Index5);
v.f = stateLGR.f(:,2);
%User Line: v = stateLGR(:,2);
u.dz0 = control.dz0; u.f = control.f;
%User Line: u = control;
%User Line: %-----------------------------------------------------------------%
%User Line: % Compute the right-hand side of the differential equations at    %
%User Line: % the N LGR points.  Each component of the right-hand side is     %
%User Line: % stored as a column vector of length N, that is each column has  %
%User Line: % the form                                                        %
%User Line: %                   [ f_i(x_1,u_1,t_1) ]                          %
%User Line: %                   [ f_i(x_2,u_2,t_2) ]                          %
%User Line: %                           .                                     %
%User Line: %                           .                                     %
%User Line: %                           .                                     %
%User Line: %                   [ f_i(x_N,u_N,t_N) ]                          %
%User Line: % where "i" is the right-hand side of the ith component of the    %
%User Line: % vector field f.  It is noted that in MATLABB the calculation of %
%User Line: % the right-hand side is vectorized.                              %
%User Line: %-----------------------------------------------------------------%
cada1td1 = zeros(80,1);
cada1td1(Gator1Data.Index6) = v.dz0;
cada1td1(Gator1Data.Index7) = u.dz0;
diffeqRHS.dz0 = cada1td1;
diffeqRHS.f = [v.f u.f];
%User Line: diffeqRHS = [v, u];
%User Line: %-----------------------------------------------------------------%
%User Line: % Compute the left-hand side of the defect constraints, recalling %
%User Line: % that the left-hand side is computed using the state at the LGR  %
%User Line: % points PLUS the final point.                                    %
%User Line: %-----------------------------------------------------------------%
cada1td1 = sparse(Gator1Data.Index8,Gator1Data.Index9,statePlusEnd.dz0,41,82);
cada1td1 = D*cada1td1;
cada1td1 = cada1td1(:);
diffeqLHS.dz0 = full(cada1td1(Gator1Data.Index10));
diffeqLHS.f = D*statePlusEnd.f;
%User Line: diffeqLHS = D*statePlusEnd;
%User Line: %-----------------------------------------------------------------%
%User Line: % Construct the defect constraints at the N LGR points.           %
%User Line: % Remember that the right-hand side needs to be scaled by the     %
%User Line: % factor (tf-t0)/2 because the rate of change of the state is     %
%User Line: % being taken with respect to $\tau\in[-1,+1]$.  Thus, we have    %
%User Line: % $dt/t\dau=(tf-t0)/2$.                                           %
%User Line: %-----------------------------------------------------------------%
cada1td1 = zeros(2,1);
cada1td1(2) = tf.dz0;
cada1td1(1) = cada1td1(1) + -t0.dz0;
cada1f1dz0 = cada1td1;
cada1f1 = tf.f - t0.f;
cada1tempdz0 = cada1f1dz0(Gator1Data.Index11);
cada1tf1 = diffeqRHS.f(Gator1Data.Index12);
cada1td1 = zeros(240,1);
cada1td1(Gator1Data.Index13) = cada1tf1(:).*cada1tempdz0;
cada1td1(Gator1Data.Index14) = cada1td1(Gator1Data.Index14) + cada1f1.*diffeqRHS.dz0;
cada1f2dz0 = cada1td1;
cada1f2 = cada1f1*diffeqRHS.f;
cada1f3dz0 = cada1f2dz0./2;
cada1f3 = cada1f2/2;
cada1td1 = zeros(640,1);
cada1td1(Gator1Data.Index15) = diffeqLHS.dz0;
cada1td1(Gator1Data.Index16) = cada1td1(Gator1Data.Index16) + -cada1f3dz0;
defects.dz0 = cada1td1;
defects.f = diffeqLHS.f - cada1f3;
%User Line: defects = diffeqLHS-(tf-t0)*diffeqRHS/2;
%User Line: %-----------------------------------------------------------------%
%User Line: % Reshape the defect contraints into a column vector.             %
%User Line: %-----------------------------------------------------------------%
cada1f1 = N.f*nstates;
defects.dz0 = defects.dz0;
defects.f = reshape(defects.f,cada1f1,1);
%User Line: defects = reshape(defects,N*nstates,1);
%User Line: %-----------------------------------------------------------------%
%User Line: % Construct the objective function plus constraint vector.        %
%User Line: %-----------------------------------------------------------------%
cada1f1dz0 = 2.*u.f(:).^(2-1).*u.dz0;
cada1f1 = u.f.^2;
cada1f2dz0 = cada1f1dz0;
cada1f2 = sum(cada1f1);
J.dz0 = 0.5.*cada1f2dz0;
J.f = 0.5*cada1f2;
%User Line: J   = (1/2)*sum(u.^(2));
cada1td1 = zeros(680,1);
cada1td1(Gator1Data.Index17) = J.dz0;
cada1td1(Gator1Data.Index18) = defects.dz0;
C.dz0 = cada1td1;
C.f = [J.f;defects.f];
%User Line: C = [J; defects];
C.dz0_size = [81,124];
C.dz0_location = Gator1Data.Index19;
end


function ADiGator_LoadData()
global ADiGator_BDFun_ADiGatorJac
ADiGator_BDFun_ADiGatorJac = load('BDFun_ADiGatorJac.mat');
return
end