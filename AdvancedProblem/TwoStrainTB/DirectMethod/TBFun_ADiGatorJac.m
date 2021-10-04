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

function C = TBFun_ADiGatorJac(z)
global ADiGator_TBFun_ADiGatorJac
if isempty(ADiGator_TBFun_ADiGatorJac); ADiGator_LoadData(); end
Gator1Data = ADiGator_TBFun_ADiGatorJac.TBFun_ADiGatorJac.Gator1Data;
% ADiGator Start Derivative Computations
%User Line: %-----------------------------------------------------------------%
%User Line: % Objective and constraint functions for the orbit-raising        %
%User Line: % problem.  This function is designed to be used with the NLP     %
%User Line: % solver SNOPT.                                                   %
%User Line: %-----------------------------------------------------------------%
%User Line: %      DO NOT FOR ANY REASON ALTER THE LINE OF CODE BELOW!        %
global psStuff nstates ncontrols npaths 
%User Line: global
global Beta1 Beta2 Beta3 Mu 
%User Line: global
global d1 d2 k1 k2 r1 r2 
%User Line: global
global p q N1 B1 B2 
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
t0.dz0 = z.dz0(327);
t0.f = z.f(t0Index.f);
%User Line: t0 = z(t0Index);
tf.dz0 = z.dz0(328);
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
S.dz0 = stateLGR.dz0(Gator1Data.Index4);
S.f = stateLGR.f(:,1);
%User Line: S = stateLGR(:,1);
T.dz0 = stateLGR.dz0(Gator1Data.Index5);
T.f = stateLGR.f(:,2);
%User Line: T = stateLGR(:,2);
L1.dz0 = stateLGR.dz0(Gator1Data.Index6);
L1.f = stateLGR.f(:,3);
%User Line: L1 = stateLGR(:,3);
L2.dz0 = stateLGR.dz0(Gator1Data.Index7);
L2.f = stateLGR.f(:,4);
%User Line: L2 = stateLGR(:,4);
I1.dz0 = stateLGR.dz0(Gator1Data.Index8);
I1.f = stateLGR.f(:,5);
%User Line: I1 = stateLGR(:,5);
I2.dz0 = stateLGR.dz0(Gator1Data.Index9);
I2.f = stateLGR.f(:,6);
%User Line: I2 = stateLGR(:,6);
u1.dz0 = control.dz0(Gator1Data.Index10);
u1.f = control.f(:,1);
%User Line: u1 = control(:,1);
u2.dz0 = control.dz0(Gator1Data.Index11);
u2.f = control.f(:,2);
%User Line: u2 = control(:,2);
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
%User Line: % Computing Derivatives
cada1f1 = Mu*N1;
cada1f2dz0 = Beta1.*S.dz0;
cada1f2 = Beta1*S.f;
cada1f3dz0 = I1.dz0./N1;
cada1f3 = I1.f/N1;
cada1td1 = zeros(80,1);
cada1td1(Gator1Data.Index12) = cada1f3(:).*cada1f2dz0;
cada1td1(Gator1Data.Index13) = cada1td1(Gator1Data.Index13) + cada1f2(:).*cada1f3dz0;
cada1f4dz0 = cada1td1;
cada1f4 = cada1f2.*cada1f3;
cada1f5dz0 = -cada1f4dz0;
cada1f5 = cada1f1 - cada1f4;
cada1f6dz0 = Beta3.*S.dz0;
cada1f6 = Beta3*S.f;
cada1f7dz0 = I2.dz0./N1;
cada1f7 = I2.f/N1;
cada1td1 = zeros(80,1);
cada1td1(Gator1Data.Index14) = cada1f7(:).*cada1f6dz0;
cada1td1(Gator1Data.Index15) = cada1td1(Gator1Data.Index15) + cada1f6(:).*cada1f7dz0;
cada1f8dz0 = cada1td1;
cada1f8 = cada1f6.*cada1f7;
cada1td1 = zeros(120,1);
cada1td1(Gator1Data.Index16) = cada1f5dz0;
cada1td1(Gator1Data.Index17) = cada1td1(Gator1Data.Index17) + -cada1f8dz0;
cada1f9dz0 = cada1td1;
cada1f9 = cada1f5 - cada1f8;
cada1f10dz0 = Mu.*S.dz0;
cada1f10 = Mu*S.f;
cada1td1 = cada1f9dz0;
cada1td1(Gator1Data.Index18) = cada1td1(Gator1Data.Index18) + -cada1f10dz0;
S_Derivative.dz0 = cada1td1;
S_Derivative.f = cada1f9 - cada1f10;
%User Line: S_Derivative = (Mu*N1) - (Beta1*S.*(I1/N1)) - (Beta3*S.*(I2/N1)) - (Mu*S);
cada1f1dz0 = r1.*u1.dz0;
cada1f1 = r1*u1.f;
cada1td1 = zeros(80,1);
cada1td1(Gator1Data.Index19) = L1.f(:).*cada1f1dz0;
cada1td1(Gator1Data.Index20) = cada1td1(Gator1Data.Index20) + cada1f1(:).*L1.dz0;
cada1f2dz0 = cada1td1;
cada1f2 = cada1f1.*L1.f;
cada1f3dz0 = Mu.*T.dz0;
cada1f3 = Mu*T.f;
cada1td1 = zeros(120,1);
cada1td1(Gator1Data.Index21) = cada1f2dz0;
cada1td1(Gator1Data.Index22) = cada1td1(Gator1Data.Index22) + -cada1f3dz0;
cada1f4dz0 = cada1td1;
cada1f4 = cada1f2 - cada1f3;
cada1f5 = p + q;
cada1f6 = r2*cada1f5;
cada1f7dz0 = -u2.dz0;
cada1f7 = 1 - u2.f;
cada1f8dz0 = -cada1f7dz0;
cada1f8 = 1 - cada1f7;
cada1f9dz0 = cada1f6.*cada1f8dz0;
cada1f9 = cada1f6*cada1f8;
cada1td1 = zeros(80,1);
cada1td1(Gator1Data.Index23) = I1.f(:).*cada1f9dz0;
cada1td1(Gator1Data.Index24) = cada1td1(Gator1Data.Index24) + cada1f9(:).*I1.dz0;
cada1f10dz0 = cada1td1;
cada1f10 = cada1f9.*I1.f;
cada1td1 = zeros(200,1);
cada1td1(Gator1Data.Index25) = cada1f4dz0;
cada1td1(Gator1Data.Index26) = cada1td1(Gator1Data.Index26) + -cada1f10dz0;
cada1f11dz0 = cada1td1;
cada1f11 = cada1f4 - cada1f10;
cada1f12dz0 = Beta2.*T.dz0;
cada1f12 = Beta2*T.f;
cada1f13dz0 = I1.dz0./N1;
cada1f13 = I1.f/N1;
cada1td1 = zeros(80,1);
cada1td1(Gator1Data.Index27) = cada1f13(:).*cada1f12dz0;
cada1td1(Gator1Data.Index28) = cada1td1(Gator1Data.Index28) + cada1f12(:).*cada1f13dz0;
cada1f14dz0 = cada1td1;
cada1f14 = cada1f12.*cada1f13;
cada1td1 = cada1f11dz0;
cada1td1(Gator1Data.Index29) = cada1td1(Gator1Data.Index29) + -cada1f14dz0;
cada1f15dz0 = cada1td1;
cada1f15 = cada1f11 - cada1f14;
cada1f16dz0 = Beta3.*T.dz0;
cada1f16 = Beta3*T.f;
cada1f17dz0 = I2.dz0./N1;
cada1f17 = I2.f/N1;
cada1td1 = zeros(80,1);
cada1td1(Gator1Data.Index30) = cada1f17(:).*cada1f16dz0;
cada1td1(Gator1Data.Index31) = cada1td1(Gator1Data.Index31) + cada1f16(:).*cada1f17dz0;
cada1f18dz0 = cada1td1;
cada1f18 = cada1f16.*cada1f17;
cada1td1 = zeros(240,1);
cada1td1(Gator1Data.Index32) = cada1f15dz0;
cada1td1(Gator1Data.Index33) = cada1td1(Gator1Data.Index33) + -cada1f18dz0;
T_Derivative.dz0 = cada1td1;
T_Derivative.f = cada1f15 - cada1f18;
%User Line: T_Derivative = (r1*u1.*L1) - (Mu*T) - (r2*(p+q)*(1-(1-u2)).*I1)                - (Beta2*T.*(I1/N1)) - (Beta3*T.*(I2/N1));
cada1f1dz0 = Beta1.*S.dz0;
cada1f1 = Beta1*S.f;
cada1f2dz0 = I1.dz0./N1;
cada1f2 = I1.f/N1;
cada1td1 = zeros(80,1);
cada1td1(Gator1Data.Index34) = cada1f2(:).*cada1f1dz0;
cada1td1(Gator1Data.Index35) = cada1td1(Gator1Data.Index35) + cada1f1(:).*cada1f2dz0;
cada1f3dz0 = cada1td1;
cada1f3 = cada1f1.*cada1f2;
cada1f4 = Mu + k1;
cada1f5dz0 = cada1f4.*L1.dz0;
cada1f5 = cada1f4*L1.f;
cada1td1 = zeros(120,1);
cada1td1(Gator1Data.Index36) = cada1f3dz0;
cada1td1(Gator1Data.Index37) = cada1td1(Gator1Data.Index37) + -cada1f5dz0;
cada1f6dz0 = cada1td1;
cada1f6 = cada1f3 - cada1f5;
cada1f7dz0 = r1.*u1.dz0;
cada1f7 = r1*u1.f;
cada1td1 = zeros(80,1);
cada1td1(Gator1Data.Index38) = L1.f(:).*cada1f7dz0;
cada1td1(Gator1Data.Index39) = cada1td1(Gator1Data.Index39) + cada1f7(:).*L1.dz0;
cada1f8dz0 = cada1td1;
cada1f8 = cada1f7.*L1.f;
cada1td1 = zeros(160,1);
cada1td1(Gator1Data.Index40) = cada1f6dz0;
cada1td1(Gator1Data.Index41) = cada1td1(Gator1Data.Index41) + -cada1f8dz0;
cada1f9dz0 = cada1td1;
cada1f9 = cada1f6 - cada1f8;
cada1f10 = p*r2;
cada1f11dz0 = -u2.dz0;
cada1f11 = 1 - u2.f;
cada1f12dz0 = cada1f10.*cada1f11dz0;
cada1f12 = cada1f10*cada1f11;
cada1td1 = zeros(80,1);
cada1td1(Gator1Data.Index42) = I1.f(:).*cada1f12dz0;
cada1td1(Gator1Data.Index43) = cada1td1(Gator1Data.Index43) + cada1f12(:).*I1.dz0;
cada1f13dz0 = cada1td1;
cada1f13 = cada1f12.*I1.f;
cada1td1 = zeros(200,1);
cada1td1(Gator1Data.Index44) = cada1f9dz0;
cada1td1(Gator1Data.Index45) = cada1td1(Gator1Data.Index45) + cada1f13dz0;
cada1f14dz0 = cada1td1;
cada1f14 = cada1f9 + cada1f13;
cada1f15dz0 = Beta2.*T.dz0;
cada1f15 = Beta2*T.f;
cada1f16dz0 = I1.dz0./N1;
cada1f16 = I1.f/N1;
cada1td1 = zeros(80,1);
cada1td1(Gator1Data.Index46) = cada1f16(:).*cada1f15dz0;
cada1td1(Gator1Data.Index47) = cada1td1(Gator1Data.Index47) + cada1f15(:).*cada1f16dz0;
cada1f17dz0 = cada1td1;
cada1f17 = cada1f15.*cada1f16;
cada1td1 = zeros(240,1);
cada1td1(Gator1Data.Index48) = cada1f14dz0;
cada1td1(Gator1Data.Index49) = cada1td1(Gator1Data.Index49) + cada1f17dz0;
cada1f18dz0 = cada1td1;
cada1f18 = cada1f14 + cada1f17;
cada1f19dz0 = Beta3.*L1.dz0;
cada1f19 = Beta3*L1.f;
cada1f20dz0 = I2.dz0./N1;
cada1f20 = I2.f/N1;
cada1td1 = zeros(80,1);
cada1td1(Gator1Data.Index50) = cada1f20(:).*cada1f19dz0;
cada1td1(Gator1Data.Index51) = cada1td1(Gator1Data.Index51) + cada1f19(:).*cada1f20dz0;
cada1f21dz0 = cada1td1;
cada1f21 = cada1f19.*cada1f20;
cada1td1 = zeros(280,1);
cada1td1(Gator1Data.Index52) = cada1f18dz0;
cada1td1(Gator1Data.Index53) = cada1td1(Gator1Data.Index53) + -cada1f21dz0;
L1_Derivative.dz0 = cada1td1;
L1_Derivative.f = cada1f18 - cada1f21;
%User Line: L1_Derivative = (Beta1*S.*(I1/N1)) - ((Mu+k1)*L1) - (r1*u1.*L1)                 + (p*r2*(1-u2).*I1) + (Beta2*T.*(I1/N1))                 - (Beta3*L1.*(I2/N1));
cada1f1 = q*r2;
cada1f2dz0 = -u2.dz0;
cada1f2 = 1 - u2.f;
cada1f3dz0 = cada1f1.*cada1f2dz0;
cada1f3 = cada1f1*cada1f2;
cada1td1 = zeros(80,1);
cada1td1(Gator1Data.Index54) = I1.f(:).*cada1f3dz0;
cada1td1(Gator1Data.Index55) = cada1td1(Gator1Data.Index55) + cada1f3(:).*I1.dz0;
cada1f4dz0 = cada1td1;
cada1f4 = cada1f3.*I1.f;
cada1f5 = Mu + k2;
cada1f6dz0 = cada1f5.*L2.dz0;
cada1f6 = cada1f5*L2.f;
cada1td1 = zeros(120,1);
cada1td1(Gator1Data.Index56) = cada1f4dz0;
cada1td1(Gator1Data.Index57) = cada1td1(Gator1Data.Index57) + -cada1f6dz0;
cada1f7dz0 = cada1td1;
cada1f7 = cada1f4 - cada1f6;
cada1td1 = zeros(80,1);
cada1td1(Gator1Data.Index58) = S.dz0;
cada1td1(Gator1Data.Index59) = cada1td1(Gator1Data.Index59) + L1.dz0;
cada1f8dz0 = cada1td1;
cada1f8 = S.f + L1.f;
cada1td1 = zeros(120,1);
cada1td1(Gator1Data.Index60) = cada1f8dz0;
cada1td1(Gator1Data.Index61) = cada1td1(Gator1Data.Index61) + T.dz0;
cada1f9dz0 = cada1td1;
cada1f9 = cada1f8 + T.f;
cada1f10dz0 = Beta3.*cada1f9dz0;
cada1f10 = Beta3*cada1f9;
cada1f11dz0 = I2.dz0./N1;
cada1f11 = I2.f/N1;
cada1tf1 = cada1f11(Gator1Data.Index62);
cada1td1 = zeros(160,1);
cada1td1(Gator1Data.Index63) = cada1tf1(:).*cada1f10dz0;
cada1td1(Gator1Data.Index64) = cada1td1(Gator1Data.Index64) + cada1f10(:).*cada1f11dz0;
cada1f12dz0 = cada1td1;
cada1f12 = cada1f10.*cada1f11;
cada1td1 = zeros(280,1);
cada1td1(Gator1Data.Index65) = cada1f7dz0;
cada1td1(Gator1Data.Index66) = cada1td1(Gator1Data.Index66) + cada1f12dz0;
L2_Derivative.dz0 = cada1td1;
L2_Derivative.f = cada1f7 + cada1f12;
%User Line: L2_Derivative = (q*r2*(1-u2).*I1) - ((Mu+k2)*L2)                 + (Beta3*(S+L1+T).*(I2/N1));
cada1f1dz0 = k1.*L1.dz0;
cada1f1 = k1*L1.f;
cada1f2 = Mu + d1;
cada1f3dz0 = cada1f2.*I1.dz0;
cada1f3 = cada1f2*I1.f;
cada1td1 = zeros(80,1);
cada1td1(Gator1Data.Index67) = cada1f1dz0;
cada1td1(Gator1Data.Index68) = cada1td1(Gator1Data.Index68) + -cada1f3dz0;
cada1f4dz0 = cada1td1;
cada1f4 = cada1f1 - cada1f3;
cada1f5dz0 = r2.*I1.dz0;
cada1f5 = r2*I1.f;
cada1td1 = cada1f4dz0;
cada1td1(Gator1Data.Index69) = cada1td1(Gator1Data.Index69) + -cada1f5dz0;
I1_Derivative.dz0 = cada1td1;
I1_Derivative.f = cada1f4 - cada1f5;
%User Line: I1_Derivative = (k1*L1) - ((Mu+d1)*I1) - (r2*I1);
cada1f1dz0 = k2.*L2.dz0;
cada1f1 = k2*L2.f;
cada1f2 = Mu + d2;
cada1f3dz0 = cada1f2.*I2.dz0;
cada1f3 = cada1f2*I2.f;
cada1td1 = zeros(80,1);
cada1td1(Gator1Data.Index70) = cada1f1dz0;
cada1td1(Gator1Data.Index71) = cada1td1(Gator1Data.Index71) + -cada1f3dz0;
I2_Derivative.dz0 = cada1td1;
I2_Derivative.f = cada1f1 - cada1f3;
%User Line: I2_Derivative = (k2*L2) - ((Mu+d2)*I2);
cada1td1 = zeros(1080,1);
cada1td1(Gator1Data.Index72) = S_Derivative.dz0;
cada1td1(Gator1Data.Index73) = T_Derivative.dz0;
cada1td1(Gator1Data.Index74) = L1_Derivative.dz0;
cada1td1(Gator1Data.Index75) = L2_Derivative.dz0;
cada1td1(Gator1Data.Index76) = I1_Derivative.dz0;
cada1td1(Gator1Data.Index77) = I2_Derivative.dz0;
diffeqRHS.dz0 = cada1td1;
diffeqRHS.f = [S_Derivative.f T_Derivative.f L1_Derivative.f L2_Derivative.f I1_Derivative.f I2_Derivative.f];
%User Line: diffeqRHS = [S_Derivative, T_Derivative, L1_Derivative,     L2_Derivative, I1_Derivative, I2_Derivative];
%User Line: %-----------------------------------------------------------------%
%User Line: % Compute the left-hand side of the defect constraints, recalling %
%User Line: % that the left-hand side is computed using the state at the LGR  %
%User Line: % points PLUS the final point.                                    %
%User Line: %-----------------------------------------------------------------%
cada1td1 = sparse(Gator1Data.Index78,Gator1Data.Index79,statePlusEnd.dz0,41,246);
cada1td1 = D*cada1td1;
cada1td1 = cada1td1(:);
diffeqLHS.dz0 = full(cada1td1(Gator1Data.Index80));
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
cada1tempdz0 = cada1f1dz0(Gator1Data.Index81);
cada1tf1 = diffeqRHS.f(Gator1Data.Index82);
cada1td1 = zeros(1560,1);
cada1td1(Gator1Data.Index83) = cada1tf1(:).*cada1tempdz0;
cada1td1(Gator1Data.Index84) = cada1td1(Gator1Data.Index84) + cada1f1.*diffeqRHS.dz0;
cada1f2dz0 = cada1td1;
cada1f2 = cada1f1*diffeqRHS.f;
cada1f3dz0 = cada1f2dz0./2;
cada1f3 = cada1f2/2;
cada1td1 = zeros(2520,1);
cada1td1(Gator1Data.Index85) = diffeqLHS.dz0;
cada1td1(Gator1Data.Index86) = cada1td1(Gator1Data.Index86) + -cada1f3dz0;
defects.dz0 = cada1td1;
defects.f = diffeqLHS.f - cada1f3;
%User Line: defects = diffeqLHS-(tf-t0)*diffeqRHS/2;
%User Line: %-----------------------------------------------------------------%
%User Line: % Construct the path constraints at the N LGR points.             %
%User Line: %-----------------------------------------------------------------%
cada1td1 = zeros(80,1);
cada1td1(Gator1Data.Index87) = S.dz0;
cada1td1(Gator1Data.Index88) = cada1td1(Gator1Data.Index88) + T.dz0;
cada1f1dz0 = cada1td1;
cada1f1 = S.f + T.f;
cada1td1 = zeros(120,1);
cada1td1(Gator1Data.Index89) = cada1f1dz0;
cada1td1(Gator1Data.Index90) = cada1td1(Gator1Data.Index90) + L1.dz0;
cada1f2dz0 = cada1td1;
cada1f2 = cada1f1 + L1.f;
cada1td1 = zeros(160,1);
cada1td1(Gator1Data.Index91) = cada1f2dz0;
cada1td1(Gator1Data.Index92) = cada1td1(Gator1Data.Index92) + L2.dz0;
cada1f3dz0 = cada1td1;
cada1f3 = cada1f2 + L2.f;
cada1td1 = zeros(200,1);
cada1td1(Gator1Data.Index93) = cada1f3dz0;
cada1td1(Gator1Data.Index94) = cada1td1(Gator1Data.Index94) + I1.dz0;
cada1f4dz0 = cada1td1;
cada1f4 = cada1f3 + I1.f;
cada1td1 = zeros(240,1);
cada1td1(Gator1Data.Index95) = cada1f4dz0;
cada1td1(Gator1Data.Index96) = cada1td1(Gator1Data.Index96) + I2.dz0;
cada1f5dz0 = cada1td1;
cada1f5 = cada1f4 + I2.f;
paths.dz0 = cada1f5dz0;
paths.f = cada1f5 - N1;
%User Line: paths = S + T + L1 + L2 + I1 + I2 - N1;
%User Line: %-----------------------------------------------------------------%
%User Line: % Reshape the defect contraints into a column vector.             %
%User Line: %-----------------------------------------------------------------%
cada1f1 = N.f*nstates;
defects.dz0 = defects.dz0;
defects.f = reshape(defects.f,cada1f1,1);
%User Line: defects = reshape(defects,N*nstates,1);
%User Line: %-----------------------------------------------------------------%
%User Line: % Reshape the Path contraints into a column vector.             %
%User Line: %-----------------------------------------------------------------%
cada1f1 = N.f*npaths;
paths.dz0 = paths.dz0;
paths.f = reshape(paths.f,cada1f1,1);
%User Line: paths = reshape(paths,N*npaths,1);
%User Line: %-----------------------------------------------------------------%
%User Line: % Construct the objective function plus constraint vector.        %
%User Line: %-----------------------------------------------------------------%
cada1td1 = zeros(80,1);
cada1td1(Gator1Data.Index97) = L2.dz0;
cada1td1(Gator1Data.Index98) = cada1td1(Gator1Data.Index98) + I2.dz0;
cada1f1dz0 = cada1td1;
cada1f1 = L2.f + I2.f;
cada1f2 = 0.5*B1;
cada1f3dz0 = 2.*u1.f(:).^(2-1).*u1.dz0;
cada1f3 = u1.f.^2;
cada1f4dz0 = cada1f2.*cada1f3dz0;
cada1f4 = cada1f2*cada1f3;
cada1td1 = zeros(120,1);
cada1td1(Gator1Data.Index99) = cada1f1dz0;
cada1td1(Gator1Data.Index100) = cada1td1(Gator1Data.Index100) + cada1f4dz0;
cada1f5dz0 = cada1td1;
cada1f5 = cada1f1 + cada1f4;
cada1f6 = 0.5*B2;
cada1f7dz0 = 2.*u2.f(:).^(2-1).*u2.dz0;
cada1f7 = u2.f.^2;
cada1f8dz0 = cada1f6.*cada1f7dz0;
cada1f8 = cada1f6*cada1f7;
cada1td1 = zeros(160,1);
cada1td1(Gator1Data.Index101) = cada1f5dz0;
cada1td1(Gator1Data.Index102) = cada1td1(Gator1Data.Index102) + cada1f8dz0;
cada1f9dz0 = cada1td1;
cada1f9 = cada1f5 + cada1f8;
J.dz0 = cada1f9dz0;
J.f = sum(cada1f9);
%User Line: J   = sum(L2 + I2 + ((1/2)*B1*u1.^2) + ((1/2)*B2*u2.^2));
cada1td1 = zeros(2920,1);
cada1td1(Gator1Data.Index103) = J.dz0;
cada1td1(Gator1Data.Index104) = defects.dz0;
cada1td1(Gator1Data.Index105) = paths.dz0;
C.dz0 = cada1td1;
C.f = [J.f;defects.f;paths.f];
%User Line: C = [J; defects; paths];
C.dz0_size = [281,328];
C.dz0_location = Gator1Data.Index106;
end


function ADiGator_LoadData()
global ADiGator_TBFun_ADiGatorJac
ADiGator_TBFun_ADiGatorJac = load('TBFun_ADiGatorJac.mat');
return
end