%-----------------------------------------------------------------%
%                     Moon Lander Problem                         %
%-----------------------------------------------------------------%
% Solve the following optimal control problem:                    %
% Minimize t_f                                                    %
% subject to the differential equation constraints                %
% dh/dt         = v                                               %
% dv/dt         = -g + u                                          %
%-----------------------------------------------------------------%
% BEGIN: DO NOT ALTER THE FOLLOWING LINES OF CODE!!!              %
%-----------------------------------------------------------------%
global psStuff nstates ncontrols g
global iGfun jGvar 
%-----------------------------------------------------------------%
% END:   DO NOT ALTER THE FOLLOWING LINES OF CODE!!!              %
%-----------------------------------------------------------------%


%-----------------------------------------------------------------%
%             Define the constants for the problem                %
%-----------------------------------------------------------------%

g = 1.5;

%-----------------------------------------------------------------%
%  Define the sizes of quantities in the optimal control problem  %
%-----------------------------------------------------------------%
nstates = 2;
ncontrols = 1;

%-----------------------------------------------------------------%
%  Define bounds on the variables in the optimal control problem  %
%-----------------------------------------------------------------%
h_0         = 10;           v_0         = -2;    
h_f         = 0;            v_f         = 0;    
hmin        = 0;            hmax        = 10;
vmin        = -5;           vmax        = 0;
umin        = 0;            umax        = 3;
t0min       = 0;            t0max       = 0;
tfmin       = 0;            tfmax       = 10000;

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
zxmin = hmin*ones(length(tau),1);
zxmax = hmax*ones(length(tau),1);
zxmin(1) = h_0; zxmax(1) = h_0;
zxmin(NLGR+1) = h_f; zxmax(NLGR+1) = h_f;

zymin = vmin*ones(length(tau),1);
zymax = vmax*ones(length(tau),1);
zymin(1) = v_0; zymax(1) = v_0;
zymin(NLGR+1) = v_f; zymax(NLGR+1) = v_f;

zumin = umin*ones(length(tau)-1,1);
zumax = umax*ones(length(tau)-1,1);

zmin = [zxmin; zymin; zumin; t0min; tfmin];
zmax = [zxmax; zymax; zumax; t0max; tfmax];

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
hguess = linspace(h_0,h_f,NLGR+1).';
vguess = linspace(v_0,v_f,NLGR+1).';
uguess = linspace(0,0,NLGR).';
t0guess = 0;
tfguess = 4;
z0 = [hguess; vguess; uguess; t0guess; tfguess];

%-----------------------------------------------------------------%
% Generate derivatives and sparsity pattern using Adigator        %
%-----------------------------------------------------------------%
% - Constraint Funtction Derivatives
xsize  = size(z0);
x      = adigatorCreateDerivInput(xsize,'z0');
output = adigatorGenJacFile('MLPFun',{x});
S_jac  = output.JacobianStructure;
[iGfun,jGvar] = find(S_jac);

% - Objective Funtcion Derivatives
xsize  = size(z0);
x      = adigatorCreateDerivInput(xsize,'z0');
output = adigatorGenJacFile('MLPObj',{x});
grd_structure = output.JacobianStructure;

%-----------------------------------------------------------------%
% Set IPOPT callback functions
%-----------------------------------------------------------------%
funcs.objective   = @(Z)MLPObj(Z);
funcs.gradient    = @(Z)MLPGrd(Z);
funcs.constraints = @(Z)MLPCon(Z);
funcs.jacobian    = @(Z)MLPJac(Z);
funcs.jacobianstructure = @()MLPJacPat(S_jac);
options.ipopt.hessian_approximation = 'limited-memory';

%-----------------------------------------------------------------%
% Set IPOPT Options %
%-----------------------------------------------------------------%
options.ipopt.tol = 1e-5;
options.ipopt.linear_solver = 'mumps';
options.ipopt.max_iter = 2000;
options.ipopt.mu_strategy = 'adaptive';
options.ipopt.ma57_automatic_scaling = 'yes';
options.ipopt.print_user_options = 'yes';
options.ipopt.output_file = ['MLP','IPOPTinfo.txt']; % print output file
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
h = z(1:NLGR+1);
v = z(NLGR+2:2*(NLGR+1));
u = z(2*(NLGR+1)+1:2*(NLGR+1)+NLGR);
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
lam_h = costate(:,1); lam_v = costate(:,2);

%-----------------------------------------------------------------%
% plot results
%-----------------------------------------------------------------%
figure(1);
plot(t,h,'-b',t,v,'-r');
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$(h(t),v(t)$','Interpreter','LaTeX');
xlabel('Time');
ylabel('States');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontName','Times','FontSize',18);
grid on;

figure(2);
plot(tLGR,u,'-o');
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$u(t)$','Interpreter','LaTeX');
xlabel('Time');
ylabel('Control');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontName','Times','FontSize',18);
grid on;

figure(4);
plot(t,lam_h,'-b',t,lam_v,'-r');
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$(\lambda_h(t),\lambda_v(t))$','Interpreter','LaTeX');
xlabel('Time');
ylabel('Co-States');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontName','Times','FontSize',18);
grid on;

