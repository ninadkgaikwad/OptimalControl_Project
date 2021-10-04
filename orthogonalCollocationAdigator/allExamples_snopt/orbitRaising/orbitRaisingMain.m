% ---------------------------------------------------%
%             Orbit-Raising Problem                  %
% ---------------------------------------------------%
% Solve the following optimal control problem:       %
% Maximize r(t_f)                                    %
% subject to the differential equation constraints   %
%   dr/dt = v_r                                      %
%   d\theta/dt = v_\theta/r                          %
%   dv_r/dt = v_\theta^2-\mu/r^2+a u_1               %
%   dv_\theta/dt = -v_rv_\theta/r+a u_2              %
% the equality path constraint                       %
%   u_1^2 + u_2^2 = 1                                %
% and the boundary conditions                        %
%   r(0) = 1                                         %
%   \theta(0) = 0                                    %
%   v_r(0) = 0                                       %
%   v_\theta(0) = 1                                  %
%   v_r(t_f) = 0                                     %
%   \sqrt{\mu/r(t_f)} = v_\theta(t_f) = 0            %
% -------------------------------------------------- %
% BEGIN: DO NOT ALTER THE FOLLOWING LINES OF CODE!!! %
% -------------------------------------------------- %
global igrid CONSTANTS psStuff nstates ncontrols npaths
% -------------------------------------------------- %
% END:   DO NOT ALTER THE FOLLOWING LINES OF CODE!!! %
% -------------------------------------------------- %

CONSTANTS.MU = 1;
CONSTANTS.m0 = 1;
CONSTANTS.mdot = 0.0749;
CONSTANTS.T = 0.1405;
nstates = 5;
ncontrols = 2;
npaths = 1;

% Bounds on State and Control
r0 = 1; theta0 = 0; vr0 = 0; vtheta0 = 1; m0 = 1;
vrf = 0;  
rmin = 0.5; rmax = 5;
thetamin = 0; thetamax = 4*pi;
vrmin = -10; vrmax = 10;
vthetamin = -10; vthetamax = 10;
mmin = 0.1; mmax = m0;
u1min = -10; u1max = 10;
u2min = -10; u2max = 10;
t0min = 0; t0max = 0;
tfmin = 3.32; tfmax = 3.32;

% Set polynomial degree in each mesh interval
N = 50;

% Set locations of mesh points
% meshPoints = linspace(-1,1,10).';
meshPoints = [-1; 1];

% Allocate Vector of Polynomial Degrees
polyDegrees = N*ones(length(meshPoints)-1,1);

% Compute LGR points, weights, and differentiation matrix.
[tau,w,D] = lgrPS(meshPoints,polyDegrees);
NLGR = length(w);
psStuff.tau = tau;
psStuff.w   = w;
psStuff.D   = D;

% Set the bounds on the NLP variables.
zrmin = rmin*ones(length(tau),1);
zrmax = rmax*ones(length(tau),1);
zrmin(1) = r0; zrmax(1) = r0;

zthetamin = thetamin*ones(length(tau),1);
zthetamax = thetamax*ones(length(tau),1);
zthetamin(1) = theta0; zthetamax(1) = theta0;

zvrmin = vrmin*ones(length(tau),1);
zvrmax = vrmax*ones(length(tau),1);
zvrmin(1) = vr0; zvrmax(1) = vr0;
zvrmin(end) = vrf; zvrmax(end) = vrf;

zvthetamin = vthetamin*ones(length(tau),1);
zvthetamax = vthetamax*ones(length(tau),1);
zvthetamin(1) = vtheta0; zvthetamax(1) = vtheta0;

zmmin = mmin*ones(length(tau),1);
zmmax = mmax*ones(length(tau),1);
zmmin(1) = m0; zmmax(1) = m0;

zu1min = u1min*ones(length(tau)-1,1);
zu1max = u1max*ones(length(tau)-1,1);

zu2min = u2min*ones(length(tau)-1,1);
zu2max = u2max*ones(length(tau)-1,1);

zmin = [zrmin; zthetamin; zvrmin; zvthetamin; zmmin; zu1min; zu2min; t0min; tfmin];
zmax = [zrmax; zthetamax; zvrmax; zvthetamax; zmmax; zu1max; zu2max; t0max; tfmax];

% Set the bounds on the NLP constraints
% There are NSTATES sets of defect constraints.
defectMin = zeros(nstates*(length(tau)-1),1);
defectMax = zeros(nstates*(length(tau)-1),1);
% There is one path constraint
pathMin = ones(length(tau)-1,1); pathMax = ones(length(tau)-1,1);
% There is one nonlinear event constraint
bcMin = 0; bcMax = 0;
objMin = -inf; objMax = inf;
Fmin = [objMin; defectMin; pathMin; bcMin];
Fmax = [objMax; defectMax; pathMax; bcMax];

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
rguess = linspace(r0,1.5,NLGR+1).';
thetaguess = linspace(theta0,theta0,NLGR+1).';
vrguess = linspace(vr0,vrf,NLGR+1).';
vthetaguess = linspace(vtheta0,vtheta0,NLGR+1).';
mguess = linspace(m0,m0,NLGR+1).';
u1guess = linspace(1,1,NLGR).';
u2guess = linspace(0,0,NLGR).';
t0guess = 0;
tfguess = 3.32;
z0 = [rguess;thetaguess;vrguess;vthetaguess;mguess;u1guess;u2guess;t0guess;tfguess];

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
userfun = 'orbitRaisingFunWrapper';

%-----------------------------------------------------------------%
% Generate derivatives and sparsity pattern using Adigator        %
%-----------------------------------------------------------------%
xsize  = size(z0);
x      = adigatorCreateDerivInput(xsize,'z0');
output = adigatorGenJacFile('orbitRaisingFun',{x});
S_jac  = output.JacobianStructure;
[iGfun,jGvar] = find(S_jac);

% Call SNOPT Using SNSOLVE Command
[z,F,info,zmul,Fmul] = snopt(z0,zmin,zmax,zmul,zstate,Fmin,Fmax,Fmul,Fstate,userfun,ObjAdd,ObjRow,AA,iAfun,jAvar,iGfun,jGvar);
% [z,F,xmul,Fmul,info,xstate,Fstate,ns,ninf,sinf,mincw,miniw,minrw] = snsolve(z0,zmin,zmax,zmul,zstate,Fmin,Fmax,Fmul,Fstate,ObjAdd,ObjRow,AA,iAfun,jAvar,iGfun,jGvar,userfun);

% Extract the state and control from the decision vector z.
% Remember that the state is approximated at the LGR points
% plus the final point, while the control is only approximated 
% at only the LGR points.
r = z(1:NLGR+1);
theta = z(NLGR+2:2*(NLGR+1));
vr = z(2*(NLGR+1)+1:3*(NLGR+1));
vtheta = z(3*(NLGR+1)+1:4*(NLGR+1));
m = z(4*(NLGR+1)+1:5*(NLGR+1));
u1 = z(5*(NLGR+1)+1:5*(NLGR+1)+NLGR);
u2 = z(5*(NLGR+1)+NLGR+1:5*(NLGR+1)+2*NLGR);
alpha = 180/pi*atan2(u1,u2);
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
lamr = costate(:,1); lamtheta = costate(:,2);
lamvr = costate(:,3); lamvtheta = costate(:,4);

figure(1);
subplot(2,2,1);
plot(t,r,'-o');
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$r(t)$','Interpreter','LaTeX');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontName','Times','FontSize',18);
grid on;

subplot(2,2,2);
plot(t,theta,'-o');
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$\theta(t)$','Interpreter','LaTeX');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontName','Times','FontSize',18);
grid on;

subplot(2,2,3);
plot(t,vr,'-o');
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$v_r(t)$','Interpreter','LaTeX');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontName','Times','FontSize',18);
grid on;

subplot(2,2,4);
plot(t,vtheta,'-o');
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$v_\theta(t)$','Interpreter','LaTeX');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontName','Times','FontSize',18);
grid on;

figure(2);
plot(tLGR,alpha,'-o');
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$\alpha(t)$','Interpreter','LaTeX');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontName','Times','FontSize',18);
grid on;

figure(3);
subplot(1,2,1);
plot(tLGR,u1,'-o');
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$u_1(t)$','Interpreter','LaTeX');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontName','Times','FontSize',18);
grid on;

subplot(1,2,2);
plot(tLGR,u2,'-o');
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$u_2(t)$','Interpreter','LaTeX');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontName','Times','FontSize',18);
grid on;

figure(4);
subplot(2,2,1);
plot(t,lamr,'-o');
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$\lambda_r(t)$','Interpreter','LaTeX');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontName','Times','FontSize',18);
grid on;

subplot(2,2,2);
plot(t,lamtheta,'-o');
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$\lambda_\theta(t)$','Interpreter','LaTeX');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontName','Times','FontSize',18);
grid on;

subplot(2,2,3);
plot(t,lamvr,'-o');
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$\lambda_{v_r}(t)$','Interpreter','LaTeX');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontName','Times','FontSize',18);
grid on;

subplot(2,2,4);
plot(t,lamvtheta,'-o');
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$\lambda_{v_\theta}(t)$','Interpreter','LaTeX');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontName','Times','FontSize',18);
grid on;


