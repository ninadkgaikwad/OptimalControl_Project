%-----------------------------------------------------------------%
%                     Ground Mobile Robot Problem                     %
%-----------------------------------------------------------------%
% Solve the following optimal control problem:                    %
% Minimize t_f                                                    %
% subject to the differential equation constraints                %
% dx/dt         = cos(Theta)                                                %
% dy/dt         = sin(Theta)                                                 %
% d(Theta)/dt   = u                                                %
%-----------------------------------------------------------------%
% BEGIN: DO NOT ALTER THE FOLLOWING LINES OF CODE!!!              %
%-----------------------------------------------------------------%
global psStuff nstates ncontrols
global iGfun jGvar
%-----------------------------------------------------------------%
% END:   DO NOT ALTER THE FOLLOWING LINES OF CODE!!!              %
%-----------------------------------------------------------------%


%-----------------------------------------------------------------%
%             Define the constants for the problem                %
%-----------------------------------------------------------------%

% No Constants for this problem

%-----------------------------------------------------------------%
%  Define the sizes of quantities in the optimal control problem  %
%-----------------------------------------------------------------%
nstates = 3;
ncontrols = 1;

%-----------------------------------------------------------------%
%  Define bounds on the variables in the optimal control problem  %
%-----------------------------------------------------------------%
x_0         = 0;    y_0         = 0;    Theta_0 = -pi;
x_f         = 0;    y_f         = 0;    Theta_f = pi;
xmin        = -100;   xmax        = 100;
ymin        = -100;   ymax        = 100;
Thetamin    = -pi;  Thetamax    = pi;
umin        = -pi;  umax        = pi;
t0min       = 0;    t0max       = 0;
tfmin       = 0;    tfmax       = 100;

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
zxmin = xmin*ones(length(tau),1);
zxmax = xmax*ones(length(tau),1);
zxmin(1) = x_0; zxmax(1) = x_0;
zxmin(NLGR+1) = x_f; zxmax(NLGR+1) = x_f;

zymin = ymin*ones(length(tau),1);
zymax = ymax*ones(length(tau),1);
zymin(1) = y_0; zymax(1) = y_0;
zymin(NLGR+1) = y_f; zymax(NLGR+1) = y_f;

zThetamin = Thetamin*ones(length(tau),1);
zThetamax = Thetamax*ones(length(tau),1);
zThetamin(1) = Theta_0; zThetamax(1) = Theta_0;

zumin = umin*ones(length(tau)-1,1);
zumax = umax*ones(length(tau)-1,1);

zmin = [zxmin; zymin; zThetamin; zumin; t0min; tfmin];
zmax = [zxmax; zymax; zThetamax; zumax; t0max; tfmax];

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
xguess = linspace(x_0,x_f,NLGR+1).';
yguess = linspace(y_0,y_f,NLGR+1).';
vguess = linspace(Theta_0,2,NLGR+1).';
uguess = linspace(0,0,NLGR).';
t0guess = 0;
tfguess = 4;
z0 = [xguess; yguess; vguess; uguess; t0guess; tfguess];

%-----------------------------------------------------------------%
% Generate derivatives and sparsity pattern using Adigator        %
%-----------------------------------------------------------------%
% - Constraint Funtction Derivatives
xsize  = size(z0);
x      = adigatorCreateDerivInput(xsize,'z0');
output = adigatorGenJacFile('GMRFun',{x});
S_jac  = output.JacobianStructure;
[iGfun,jGvar] = find(S_jac);

% - Objective Funtcion Derivatives
xsize  = size(z0);
x      = adigatorCreateDerivInput(xsize,'z0');
output = adigatorGenJacFile('GMRObj',{x});
grd_structure = output.JacobianStructure;

%-----------------------------------------------------------------%
% Set IPOPT callback functions
%-----------------------------------------------------------------%
funcs.objective   = @(Z)GMRObj(Z);
funcs.gradient    = @(Z)GMRGrd(Z);
funcs.constraints = @(Z)GMRCon(Z);
funcs.jacobian    = @(Z)GMRJac(Z);
funcs.jacobianstructure = @()GMRJacPat(S_jac);
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
options.ipopt.output_file = ['brachistochrone','IPOPTinfo.txt']; % print output file
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
x = z(1:NLGR+1);
y = z(NLGR+2:2*(NLGR+1));
Theta = z(2*(NLGR+1)+1:3*(NLGR+1));
u = z(3*(NLGR+1)+1:3*(NLGR+1)+NLGR);
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
lam_x = costate(:,1); lam_y = costate(:,2); lam_Theta = costate(:,3);

%-----------------------------------------------------------------%
% plot results
%-----------------------------------------------------------------%
figure(1);
plot(t,x,'-bs',t,y,'-ro',t,Theta','-gd');
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$(x(t),y(t),\theta(t))$','Interpreter','LaTeX');
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

figure(3);
plot(x,y,'-o');
xl = xlabel('$x(t)$','Interpreter','LaTeX');
yl = ylabel('$y(t)$','Interpreter','LaTeX');
xlabel('Time');
ylabel('States');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontName','Times','FontSize',18);
grid on;

figure(4);
plot(t,lam_x,'-bs',t,lam_y,'-ro',t,lam_Theta,'-gd');
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$(\lambda_x(t),\lambda_y(t),\lambda_{\theta}(t))$','Interpreter','LaTeX');
xlabel('Time');
ylabel('Co-States');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontName','Times','FontSize',18);
grid on;

