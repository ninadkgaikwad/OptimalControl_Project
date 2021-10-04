%-----------------------------------------------------------------%
%               Two Strain Tuberculosis Strain Problem            %
%-----------------------------------------------------------------%
% Solve the following optimal control problem:                    %
% Minimize t_f                                                    %
% subject to the differential equation constraints                %
% dx/dt         = v                                               %
% dv/dt         = u                                          %

clear all;
close all;
clc;

%-----------------------------------------------------------------%
% BEGIN: DO NOT ALTER THE FOLLOWING LINES OF CODE!!!              %
%-----------------------------------------------------------------%
global psStuff nstates ncontrols npaths 
global iGfun jGvar 
global Beta1 Beta2 Beta3 Mu
global d1 d2 k1 k2 r1 r2 
global p q N1 B1 B2 
%-----------------------------------------------------------------%
% END:   DO NOT ALTER THE FOLLOWING LINES OF CODE!!!              %
%-----------------------------------------------------------------%


%-----------------------------------------------------------------%
%             Define the constants for the problem                %
%-----------------------------------------------------------------%

Beta1 = 13;
Beta2 = 13;
Beta3 = 0.029;
Mu = 0.0143;
d1 = 0;
d2 = 0;
k1 = 0.5;
k2 = 1;
r1 = 2;
r2 = 1;
p = 0.4;
q = 0.1;
N1 = 30000;
B1 = 50;
B2 = 500;

%-----------------------------------------------------------------%
%  Define the sizes of quantities in the optimal control problem  %
%-----------------------------------------------------------------%
nstates = 6;
ncontrols = 2;
npaths = 1;

%-----------------------------------------------------------------%
%  Define bounds on the variables in the optimal control problem  %
%-----------------------------------------------------------------%

% Know Initial States
S_0             = (76*N1)/(120);
T_0             = (N1)/(120); 
L1_0            = (36*N1)/(120); 
L2_0            = (2*N1)/(120); 
I1_0            = (4*N1)/(120); 
I2_0            = (N1)/(120); 

% Min and Max for States
Smin         = 0;               Smax         = 30000;
Tmin         = 0;               Tmax         = 30000;
L1min        = 0;               L1max        = 30000;
L2min        = 0;               L2max        = 30000;
I1min        = 0;               I1max        = 30000;
I2min        = 0;               I2max        = 30000;

% Min and Max for Controls
u1min        = 0.05;            u1max        = 0.95;
u2min        = 0.05;            u2max        = 0.95;

% Min and Max for Time
t0min       = 0;             t0max       = 0;
tfmin       = 5;             tfmax       = 5;

%-----------------------------------------------------------------%
% In this section, we define the three type of discretizations    %
% that can be employed.  These three approaches are as follows:   %
%    (1) p-method  = global pseudospectral method                 %
%                  = single interval and the degree of the        %
%                    polynomial in the interval can be varied     %
%    (2) h-method  = fixed-degree polynomial in each interval     %
%                    and the number of intervals can be varied    %
%    (3) hp-method = can vary BOTH the degree of the polynomial   %
%                    in each interval and the number of intervals %
%                                                                 %
% For simplicity in this tutorial, we will allow for either a     %
% p-method or an h-method.  Regardless of which method is being   %
% employed, the user needs to specify the following parameters:   %
%    (a) N = Polynomial Degree                                    %
%    (b) meshPoints = Set of Monotonically Increasing Mesh Points %
%                     on the Interval $\tau\in[-1,+1]$.           %
%                                                                 %
% When using a p-method, the parameters N and meshPoints must be  %
% specified as follows:                                           %
%    (i)  meshPoints = [-1 1]                                     %
%    (ii) N = Choice of Polynomial Degree (e.g., N=10, N=20)      %
% When using an h-method, the parameters N and meshPoints must be %
% specified as follows:                                           %
%    (i)  meshPoints = $[\tau_1,\tau_2,\tau_3,\ldots,\tau_N]$     %
%                      where $\tau_1 = -1$, $\tau_N = 1$ and      %
%                      (\tau_2,\ldots,\tau_{N-1}) are             %
%                      monotonically increasing on the open       %
%                      interval $(-1,+1)$.                        %
%-----------------------------------------------------------------%
%      Compute Points, Weights, and Differentiation Matrix        %
%-----------------------------------------------------------------%
%-----------------------------------------------------------------%
% Choose Polynomial Degree and Number of Mesh Intervals           %
% numIntervals = 1 ===> p-method                                  %
% numIntervals > 1 ===> h-method                                  %
%-----------------------------------------------------------------%
N = 4;
numIntervals = 10;
%-----------------------------------------------------------------%
% DO NOT ALTER THE LINE OF CODE SHOWN BELOW!                      %
%-----------------------------------------------------------------%
meshPoints = linspace(-1,1,numIntervals+1).';  
polyDegrees = N*ones(numIntervals,1);
[tau,w,D] = lgrPS(meshPoints,polyDegrees);
psStuff.tau = tau; psStuff.w = w; psStuff.D = D; NLGR = length(w);
%-----------------------------------------------------------------%
% DO NOT ALTER THE LINES OF CODE SHOWN ABOVE!                     %
%-----------------------------------------------------------------%

%-----------------------------------------------------------------%
% Set the bounds on the variables in the NLP.                     %
%-----------------------------------------------------------------%

% Min and Max for States
zSmin = Smin*ones(length(tau),1);
zSmax = Smax*ones(length(tau),1);
zSmin(1) = S_0; zSmax(1) = S_0;

zTmin = Tmin*ones(length(tau),1);
zTmax = Tmax*ones(length(tau),1);
zTmin(1) = T_0; zTmax(1) = T_0;

zL1min = L1min*ones(length(tau),1);
zL1max = L1max*ones(length(tau),1);
zL1min(1) = L1_0; zL1max(1) = L1_0;

zL2min = L2min*ones(length(tau),1);
zL2max = L2max*ones(length(tau),1);
zL2min(1) = L2_0; zL2max(1) = L2_0;

zI1min = I1min*ones(length(tau),1);
zI1max = I1max*ones(length(tau),1);
zI1min(1) = I1_0; zI1max(1) = I1_0;

zI2min = I2min*ones(length(tau),1);
zI2max = I2max*ones(length(tau),1);
zI2min(1) = I2_0; zI2max(1) = I2_0;

% Min and Max for Controls
zu1min = u1min*ones(length(tau)-1,1);
zu1max = u1max*ones(length(tau)-1,1);
zu2min = u2min*ones(length(tau)-1,1);
zu2max = u2max*ones(length(tau)-1,1);

zmin = [zSmin; zTmin; zL1min; zL2min; zI1min; zI2min; zu1min; ...
    zu2min; t0min; tfmin];
zmax = [zSmax; zTmax; zL1max; zL2max; zI1max; zI2max; zu1max; ...
    zu2max; t0max; tfmax];

%-----------------------------------------------------------------%
% Set the bounds on the constraints in the NLP.                   %
%-----------------------------------------------------------------%
defectMin = zeros(nstates*(length(tau)-1),1);
defectMax = zeros(nstates*(length(tau)-1),1);
pathMin = zeros(length(tau)-1,1); pathMax = zeros(length(tau)-1,1);
eventMin = []; eventMax = [];
objMin = 0; objMax = inf;
Fmin = [objMin; defectMin; pathMin; eventMin];
Fmax = [objMax; defectMax; pathMax; eventMax];

%-----------------------------------------------------------------%
% Supply an initial guess for the NLP.                            %
%-----------------------------------------------------------------%
% Sguess = linspace(Smin,Smax,NLGR+1).';
% Tguess = linspace(Tmin,Tmax,NLGR+1).';
% L1guess = linspace(L1min,L1max,NLGR+1).';
% L2guess = linspace(L2min,L2max,NLGR+1).';
% I1guess = linspace(I1min,I1max,NLGR+1).';
% I2guess = linspace(I2min,I2max,NLGR+1).';
% u1guess = linspace(u1min,u1max,NLGR).';
% u2guess = linspace(u2min,u2max,NLGR).';
Sguess = S_0*ones(NLGR+1,1);
Tguess = T_0*ones(NLGR+1,1);
L1guess = L1_0*ones(NLGR+1,1);
L2guess = L2_0*ones(NLGR+1,1);
I1guess = I1_0*ones(NLGR+1,1);
I2guess = I2_0*ones(NLGR+1,1);
u1guess = ((u1min+u1max)/2)*ones(NLGR,1);
u2guess = ((u2min+u2max)/2)*ones(NLGR,1);

t0guess = 0;
tfguess = 5;
z0 = [Sguess; Tguess;  L1guess; L2guess; I1guess; I2guess; u1guess;...
    u2guess; t0guess; tfguess];

%-----------------------------------------------------------------%
% Generate derivatives and sparsity pattern using Adigator        %
%-----------------------------------------------------------------%
% - Constraint Funtction Derivatives
xsize  = size(z0);
x      = adigatorCreateDerivInput(xsize,'z0');
output = adigatorGenJacFile('TBFun',{x});
S_jac  = output.JacobianStructure;
[iGfun,jGvar] = find(S_jac);

% - Objective Funtcion Derivatives
xsize  = size(z0);
x      = adigatorCreateDerivInput(xsize,'z0');
output = adigatorGenJacFile('TBObj',{x});
grd_structure = output.JacobianStructure;

%-----------------------------------------------------------------%
% Set IPOPT callback functions
%-----------------------------------------------------------------%
funcs.objective   = @(Z)TBObj(Z);
funcs.gradient    = @(Z)TBGrd(Z);
funcs.constraints = @(Z)TBCon(Z);
funcs.jacobian    = @(Z)TBJac(Z);
funcs.jacobianstructure = @()TBJacPat(S_jac);
options.ipopt.hessian_approximation = 'limited-memory';

%-----------------------------------------------------------------%
% Set IPOPT Options %
%-----------------------------------------------------------------%
options.ipopt.tol = 1e-5;
options.ipopt.linear_solver = 'ma57';
options.ipopt.max_iter = 10000;
options.ipopt.mu_strategy = 'adaptive';
options.ipopt.ma57_automatic_scaling = 'yes';
options.ipopt.print_user_options = 'yes';
options.ipopt.output_file = ['TB','IPOPTinfo.txt']; % print output file
options.ipopt.print_level = 5; % set print level default

options.lb = zmin; % Lower bound on the variables.
options.ub = zmax; % Upper bound on the variables.
options.cl = Fmin; % Lower bounds on the constraint functions.
options.cu = Fmax; % Upper bounds on the constraint functions.

%-----------------------------------------------------------------%
% Call IPOPT
%-----------------------------------------------------------------%
[z, info] = ipopt(z0,funcs,options);

%-----------------------------------------------------------------%
% extract lagrange multipliers from ipopt output, info
%-----------------------------------------------------------------%
Fmul = info.lambda;

%-----------------------------------------------------------------%
% Extract the state and control from the decision vector z.       %
% Remember that the state is approximated at the LGR points       %
% plus the final point, while the control is only approximated    %
% at only the LGR points.                                         %
%-----------------------------------------------------------------%
S = z(1:NLGR+1);
T = z(NLGR+2:2*(NLGR+1));
L1 = z(2*(NLGR+1)+1:3*(NLGR+1));
L2 = z(3*(NLGR+1)+1:4*(NLGR+1));
I1 = z(4*(NLGR+1)+1:5*(NLGR+1));
I2 = z(5*(NLGR+1)+1:6*(NLGR+1));
u1 = z(6*(NLGR+1)+1:6*(NLGR+1)+NLGR);
u2 = z(6*(NLGR+1)+NLGR+1:6*(NLGR+1)+2*NLGR);
t0 = z(end-1);
tf = z(end);
t = (tf-t0)*(tau+1)/2+t0; % Time for Plotting States
tLGR = t(1:end-1); % Time for Plotting Control

%-----------------------------------------------------------------%
% Extract the Lagrange multipliers corresponding                  %
% the defect constraints.                                         %
%-----------------------------------------------------------------%
multipliersDefects = Fmul(2:nstates*NLGR+1);
multipliersDefects = reshape(multipliersDefects,NLGR,nstates);
%-----------------------------------------------------------------%
% Compute the costates at the LGR points via transformation       %
%-----------------------------------------------------------------%
costateLGR = inv(diag(w))*multipliersDefects;
%-----------------------------------------------------------------%
% Compute the costate at the tau=+1 via transformation            %
%-----------------------------------------------------------------%
costateF = D(:,end).'*multipliersDefects;
%-----------------------------------------------------------------%
% Now assemble the costates into a single matrix                  %
%-----------------------------------------------------------------%
costate = [costateLGR; costateF];
lam_S = costate(:,1); lam_T = costate(:,2);
lam_L1 = costate(:,3); lam_L2 = costate(:,4);
lam_I1 = costate(:,5); lam_I2 = costate(:,6);

%-----------------------------------------------------------------%
% plot results
%-----------------------------------------------------------------%
figure(1);

subplot(2,3,1)
plot(t,S,'-b');
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$(S(t)$','Interpreter','LaTeX');
xlabel('Time');
ylabel('States - S');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontName','Times','FontSize',18);
grid on;

%figure(2);
subplot(2,3,2)
plot(t,T,'-b');
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$T(t)$','Interpreter','LaTeX');
xlabel('Time');
ylabel('State - T');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontName','Times','FontSize',18);
grid on;

%figure(3);
subplot(2,3,3)
plot(t,L1,'-b');
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$L_{1}(t)$','Interpreter','LaTeX');
xlabel('Time');
ylabel('State - L1');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontName','Times','FontSize',18);
grid on;

%figure(4);
subplot(2,3,4)
plot(t,L2,'-b');
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$L_{2}(t)$','Interpreter','LaTeX');
xlabel('Time');
ylabel('State - L2');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontName','Times','FontSize',18);
grid on;

%figure(5);
subplot(2,3,5)
plot(t,I1,'-b');
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$I_{1}(t)$','Interpreter','LaTeX');
xlabel('Time');
ylabel('State - I1');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontName','Times','FontSize',18);
grid on;

%figure(6);
subplot(2,3,6)
plot(t,I2,'-b');
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$I_{2}(t)$','Interpreter','LaTeX');
xlabel('Time');
ylabel('State - I2');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontName','Times','FontSize',18);
grid on;

figure(2);

%figure(7);
subplot(1,2,1)
plot(tLGR,u1,'-o');
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$u_{1}(t)$','Interpreter','LaTeX');
xlabel('Time');
ylabel('Control - U1');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontName','Times','FontSize',18);
grid on;

%figure(8);
subplot(1,2,2)
plot(tLGR,u2,'-o');
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$u_{2}(t)$','Interpreter','LaTeX');
xlabel('Time');
ylabel('Control - U2');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontName','Times','FontSize',18);
grid on;

figure(3)
%figure(9);
plot(t,lam_S,'-b',t,lam_T,'-r',t,lam_L1,'-g',t,lam_L2,'-k'...
    ,t,lam_I1,'-c',t,lam_I2,'-m');
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$(\lambda_{S}(t),\lambda_{T}(t),\lambda_{L_1}(t),\lambda_{L_2}(t),\lambda_T{I_1}t),\lambda_{I_2}(t))$',...
    'Interpreter','LaTeX');
xlabel('Time');
ylabel('Co-States');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontName','Times','FontSize',18);
grid on;

