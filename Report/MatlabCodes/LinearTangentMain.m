%-----------------------------------------------------------------%
%                     Linear Tangent Steering Problem             %

clear all;
close all;
clc;

%-----------------------------------------------------------------%
% Solve the following optimal control problem:                    %
% Minimize t_f                                                    %
% subject to the differential equation constraints                %
% dx1/dt         = x3                                             %
% dx2/dt         = x4                                             %
% dx3/dt         = a*cos(u)                                       %
% dx4/dt         = a*sin(u)                                       %
%-----------------------------------------------------------------%
% BEGIN: DO NOT ALTER THE FOLLOWING LINES OF CODE!!!              %
%-----------------------------------------------------------------%
global psStuff nstates ncontrols a
global iGfun jGvar 
%-----------------------------------------------------------------%
% END:   DO NOT ALTER THE FOLLOWING LINES OF CODE!!!              %
%-----------------------------------------------------------------%


%-----------------------------------------------------------------%
%             Define the constants for the problem                %
%-----------------------------------------------------------------%

a = 100;

%-----------------------------------------------------------------%
%  Define the sizes of quantities in the optimal control problem  %
%-----------------------------------------------------------------%
nstates = 4;
ncontrols = 1;

%-----------------------------------------------------------------%
%  Define bounds on the variables in the optimal control problem  %
%-----------------------------------------------------------------%
x1_0         = 0;              

x2_0         = 0;               
x2_f         = 5; 

x3_0         = 0;               
x3_f         = 45; 

x4_0         = 0;               
x4_f         = 0; 

x1min        = 0;             x1max        = 100;
x2min        = 0;             x2max        = 100;
x3min        = 0;             x3max        = 100;
x4min        = 0;             x4max        = 100;

umin         = -pi/2;         umax         = pi/2;

t0min       = 0;              t0max       = 0;
tfmin       = 0;              tfmax       = 100;

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
numIntervals = 30;
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
zx1min = x1min*ones(length(tau),1);
zx1max = x1max*ones(length(tau),1);
zx1min(1) = x1_0; zx1max(1) = x1_0;

zx2min = x2min*ones(length(tau),1);
zx2max = x2max*ones(length(tau),1);
zx2min(1) = x2_0; zx2max(1) = x2_0;
zx2min(NLGR+1) = x2_f; zx2max(NLGR+1) = x2_f;

zx3min = x3min*ones(length(tau),1);
zx3max = x3max*ones(length(tau),1);
zx3min(1) = x3_0; zx3max(1) = x3_0;
zx3min(NLGR+1) = x3_f; zx3max(NLGR+1) = x3_f;

zx4min = x4min*ones(length(tau),1);
zx4max = x4max*ones(length(tau),1);
zx4min(1) = x4_0; zx4max(1) = x4_0;
zx4min(NLGR+1) = x4_f; zx4max(NLGR+1) = x4_f;

zumin = umin*ones(length(tau)-1,1);
zumax = umax*ones(length(tau)-1,1);

zmin = [zx1min; zx2min; zx3min; zx4min; zumin; t0min; tfmin];
zmax = [zx1max; zx2max; zx3max; zx4max; zumax; t0max; tfmax];

%-----------------------------------------------------------------%
% Set the bounds on the constraints in the NLP.                   %
%-----------------------------------------------------------------%
defectMin = zeros(nstates*(length(tau)-1),1);
defectMax = zeros(nstates*(length(tau)-1),1);
pathMin = []; pathMax = [];
eventMin = []; eventMax = [];
objMin = 0; objMax = inf;
Fmin = [objMin; defectMin; pathMin; eventMin];
Fmax = [objMax; defectMax; pathMax; eventMax];

%-----------------------------------------------------------------%
% Supply an initial guess for the NLP.                            %
%-----------------------------------------------------------------%
x1guess = x1_0*ones(NLGR+1,1)+randn(NLGR+1,1);
x2guess = x2_0*ones(NLGR+1,1)+randn(NLGR+1,1);
x3guess = x3_0*ones(NLGR+1,1)+randn(NLGR+1,1);
x4guess = x4_0*ones(NLGR+1,1)+randn(NLGR+1,1);
uguess = ((umin-umax)/2)*ones(NLGR,1)+randn(NLGR,1);
t0guess = 0;
tfguess = 1+randn(1,1);

z0 = [x1guess; x2guess; x3guess; x4guess; uguess; t0guess; tfguess];

%-----------------------------------------------------------------%
% Generate derivatives and sparsity pattern using Adigator        %
%-----------------------------------------------------------------%
% - Constraint Funtction Derivatives
xsize  = size(z0);
x      = adigatorCreateDerivInput(xsize,'z0');
output = adigatorGenJacFile('LinearTangentFun',{x});
S_jac  = output.JacobianStructure;
[iGfun,jGvar] = find(S_jac);

% - Objective Funtcion Derivatives
xsize  = size(z0);
x      = adigatorCreateDerivInput(xsize,'z0');
output = adigatorGenJacFile('LinearTangentObj',{x});
grd_structure = output.JacobianStructure;

%-----------------------------------------------------------------%
% Set IPOPT callback functions
%-----------------------------------------------------------------%
funcs.objective   = @(Z)LinearTangentObj(Z);
funcs.gradient    = @(Z)LinearTangentGrd(Z);
funcs.constraints = @(Z)LinearTangentCon(Z);
funcs.jacobian    = @(Z)LinearTangentJac(Z);
funcs.jacobianstructure = @()LinearTangentJacPat(S_jac);
options.ipopt.hessian_approximation = 'limited-memory';

%-----------------------------------------------------------------%
% Set IPOPT Options %
%-----------------------------------------------------------------%
options.ipopt.tol = 1e-5;
options.ipopt.linear_solver = 'ma57';
options.ipopt.max_iter = 20000;
options.ipopt.mu_strategy = 'adaptive';
options.ipopt.ma57_automatic_scaling = 'yes';
options.ipopt.print_user_options = 'yes';
options.ipopt.output_file = ['LinearTangent','IPOPTinfo.txt']; % print output file
options.ipopt.print_level = 5; % set print level default

options.lb = zmin; % Lower bound on the variables.
options.ub = zmax; % Upper bound on the variables.
options.cl = Fmin; % Lower bounds on the constraint functions.
options.cu = Fmax; % Upper bounds on the constraint functions.

%-----------------------------------------------------------------%
% Call IPOPT
%-----------------------------------------------------------------%
tic; % Timing the Process

[z, info] = ipopt(z0,funcs,options);

TimeTaken = toc;

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
x1 = z(1:NLGR+1);
x2 = z(NLGR+2:2*(NLGR+1));
x3 = z(2*(NLGR+1)+1:3*(NLGR+1));
x4 = z(3*(NLGR+1)+1:4*(NLGR+1));
u = z(4*(NLGR+1)+1:4*(NLGR+1)+NLGR);
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
lam_x1 = costate(:,1); lam_x2 = costate(:,2);
lam_x3 = costate(:,3); lam_x4 = costate(:,4);

%-----------------------------------------------------------------%
% plot results
%-----------------------------------------------------------------%
% Plotting States 
figure(1)
hold on
grid on
plot(t,x1,'-g','LineWidth',1.5);
plot(t,x2,'-.m','LineWidth',1.5);
plot(t,x3,'-r','LineWidth',1.5);
plot(t,x4,'--b','LineWidth',1.5);
title('States vs. Time','Interpreter','latex');
xlabel('Time (sec)','Interpreter','latex');
ylabel('States','Interpreter','latex');
legend1=legend('$x_{1}$','$x_{2}$','$x_{3}$','$x_{4}$');
set(legend1,'Interpreter','latex');
hold off;

% Plotting Control
figure(2)
hold on
grid on
plot(tLGR,u,'-k','LineWidth',1.5);
title('Control - $u(t)$ vs. Time','Interpreter','latex');
xlabel('Time (sec)','Interpreter','latex');
ylabel('$u(t)$','Interpreter','latex');
hold off;

% Plotting Costates
figure(3)
hold on
grid on
plot(t,lam_x1,'-g','LineWidth',1.5);
plot(t,lam_x2,'-.m','LineWidth',1.5);
plot(t,lam_x3,'-r','LineWidth',1.5);
plot(t,lam_x4,'--b','LineWidth',1.5);
title('Co-States vs. Time','Interpreter','latex');
xlabel('Time (sec)','Interpreter','latex');
ylabel('Co-States','Interpreter','latex');
legend2=legend('$\lambda_{x_{1}}$','$\lambda_{x_{2}}$',...
    '$\lambda_{x_{3}}$','$\lambda_{x_{4}}$');
set(legend2,'Interpreter','latex');
hold off;

fprintf('Time taken to solve the Collocation Problem = %.4f',TimeTaken)



