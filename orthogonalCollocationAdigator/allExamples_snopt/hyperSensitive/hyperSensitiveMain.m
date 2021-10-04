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
global igrid gravity psStuff nstates ncontrols
global iGfun jGvar
% -------------------------------------------------- %
% END:   DO NOT ALTER THE FOLLOWING LINES OF CODE!!! %
% -------------------------------------------------- %

% Provide value for gravity
nstates = 1;
ncontrols = 1;

% Bounds on State and Control
x0 = 1; xf = 1; tf = 50;
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
N = 4;
numIntervals = 40;
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

% Set the derivative option
snseti('Derivative Option',1);
% Set the derivative verification level
snseti('Verify Level',-1);
% Set name of SNOPT print file
snprint('snoptmain.out');
% Print CPU times at bottom of SNOPT print file
snseti('Timing level',3);
% Echo SNOPT Output to MATLAB Command Window 
snscreen on;

% Supply an initial guess
xguess = linspace(x0,xf,NLGR+1).';
uguess = linspace(0,0,NLGR).';
t0guess = 0;
tfguess = 50;
z0 = [xguess; uguess; t0guess; tfguess];

% Set initial guess on basis and Lagrange multipliers to zero
zmul = zeros(size(z0));
zstate = zmul;
Fmul = zeros(size(Fmin));
Fstate = Fmul;
ObjAdd = 0;
% Row 1 is the objective function row
ObjRow = 1;
% Assume for simplicity that all constraints 
% are purely nonlinear
AA = [];
iAfun = [];
jAvar = [];
userfun = 'hyperSensitiveFunWrapper';

%-----------------------------------------------------------------%
% Generate derivatives and sparsity pattern using Adigator        %
%-----------------------------------------------------------------%
xsize  = size(z0);
x      = adigatorCreateDerivInput(xsize,'z0');
output = adigatorGenJacFile('hyperSensitiveFun',{x});
S_jac  = output.JacobianStructure;
[iGfun,jGvar] = find(S_jac);


% Call SNOPT Using SNSOLVE Command
[z,F,info,zmul,Fmul] = snopt(z0,zmin,zmax,zmul,zstate,Fmin,Fmax,Fmul,Fstate,userfun,ObjAdd,ObjRow,AA,iAfun,jAvar,iGfun,jGvar);
% [z,F,xmul,Fmul,info,xstate,Fstate,ns,ninf,sinf,mincw,miniw,minrw] = snsolve(z0,zmin,zmax,zmul,zstate,Fmin,Fmax,Fmul,Fstate,ObjAdd,ObjRow,AA,iAfun,jAvar,iGfun,jGvar,userfun);

% Extract the state and control from the decision vector z.
% Remember that the state is approximated at the LGR points
% plus the final point, while the control is only approximated 
% at only the LGR points.
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

