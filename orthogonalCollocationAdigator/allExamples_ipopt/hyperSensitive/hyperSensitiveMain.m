% ---------------------------------------------------%
%            Brachistochrone Problem                 %
% ---------------------------------------------------%
% Solve the following optimal control problem:       %
% Minimize                                           %
%   0.5\int_{0}^{t_f} (x^2+u^2)dt                    %
% subject to the differential equation constraint    %
%   dx/dt = -x^3 + u                                 %
% and the boundary conditions                        %
%   x(0) = x_0, x(t_f) = x_f                         %
% ---------------------------------------------------%

% -------------------------------------------------- %
% BEGIN: DO NOT ALTER THE FOLLOWING LINES OF CODE!!! %
% -------------------------------------------------- %
global igrid psStuff nstates ncontrols
global iGfun jGvar
% -------------------------------------------------- %
% END:   DO NOT ALTER THE FOLLOWING LINES OF CODE!!! %
% -------------------------------------------------- %

% Provide value for gravity
nstates = 1;
ncontrols = 1;

% Bounds on State and Control
x0 = 1.5; xf = 1; tf = 50;
xmin = -50; xmax = 50;
umin = -50; umax = 50;
t0min = 0; t0max = 0;
tfmin = tf; tfmax = tf;

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
N = 10;
numIntervals = 50;
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

% Set the bounds on the NLP variables.
zxmin = xmin*ones(length(tau),1);
zxmax = xmax*ones(length(tau),1);
zxmin(1) = x0; zxmax(1) = x0;
zxmin(NLGR+1) = xf; zxmax(NLGR+1) = xf;
zumin = umin*ones(length(tau)-1,1);
zumax = umax*ones(length(tau)-1,1);

zmin = [zxmin; zumin; t0min; tfmin];
zmax = [zxmax; zumax; t0max; tfmax];

% Set the bounds on the NLP constraints
defectMin = zeros(nstates*(length(tau)-1),1);
defectMax = zeros(nstates*(length(tau)-1),1);
pathMin = []; pathMax = [];
eventMin = []; eventMax = [];
objMin = 0; objMax = inf;
Fmin = [objMin; defectMin; pathMin; eventMin];
Fmax = [objMax; defectMax; pathMax; eventMax];

% Supply an initial guess
xguess = linspace(x0,xf,NLGR+1).';
uguess = linspace(0,0,NLGR).';
t0guess = 0;
tfguess = tf;
z0 = [xguess; uguess; t0guess; tfguess];

%-----------------------------------------------------------------%
% Generate derivatives and sparsity pattern using Adigator        %
%-----------------------------------------------------------------%
% - Constraint function derivatives
xsize  = size(z0);
x      = adigatorCreateDerivInput(xsize,'z0');
output = adigatorGenJacFile('hyperSensitiveFun',{x});
S_jac  = output.JacobianStructure;
[iGfun,jGvar] = find(S_jac);

% - Objective function derivatives
xsize  = size(z0);
x      = adigatorCreateDerivInput(xsize,'z0');
output = adigatorGenJacFile('hyperSensitiveObj',{x});
grd_structure = output.JacobianStructure;

%-----------------------------------------------------------------%
% set IPOPT callback functions
%-----------------------------------------------------------------%
funcs.objective   = @(Z)hyperSensitiveObj(Z);
funcs.gradient    = @(Z)hyperSensitiveGrd(Z);
funcs.constraints = @(Z)hyperSensitiveCon(Z);
funcs.jacobian    = @(Z)hyperSensitiveJac(Z);
funcs.jacobianstructure = @()hyperSensitiveJacPat(S_jac);
options.ipopt.hessian_approximation = 'limited-memory';

%-----------------------------------------------------------------%
% Set IPOPT Options %
%-----------------------------------------------------------------%
options.ipopt.tol = 1e-8;
options.ipopt.linear_solver = 'ma57';
options.ipopt.max_iter = 2000;
options.ipopt.mu_strategy = 'adaptive';
options.ipopt.ma57_automatic_scaling = 'yes';
options.ipopt.print_user_options = 'yes';
options.ipopt.output_file = ['hyperSensitive','IPOPTinfo.txt']; % print output file
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
% Extract the state and control from the decision vector z.
% Remember that the state is approximated at the LGR points
% plus the final point, while the control is only approximated 
% at only the LGR points.
%-----------------------------------------------------------------%
x = z(1:NLGR+1);
u = z((NLGR+1)+1:(NLGR+1)+NLGR);
t0 = z(end-1);
tf = z(end);
t = (tf-t0)*(tau+1)/2+t0;
tLGR = t(1:end-1);

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
lamx = costate;

%--------------%
% Plot Results %
%--------------%
figure(1)
plot(t,x,'k-o',t,lamx,'r-o');
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$x(t),\lambda_x(t)$','Interpreter','LaTeX');
ll = legend('$x(t)$','$\lambda_x(t)$','Interpreter','LaTeX');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(ll,'FontSize',18);
set(gca,'FontName','Times','FontSize',18);
grid on;

figure(2)
plot(x,lamx,'-o');
xl = xlabel('$x$','Interpreter','LaTeX');
yl = ylabel('$\lambda_x$','Interpreter','LaTeX');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontName','Times','FontSize',18);
grid on;

figure(3)
plot(tLGR,u,'-o');
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$u(t)$','Interpreter','LaTeX');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontName','Times','FontSize',18);
grid on;

