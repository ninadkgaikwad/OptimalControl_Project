function [ Error ] = LinearTangent_ODES0lver( X0_RootFinder )

%% ODE Solver and Error Calculator : Problem 1 - Part 2

%% Inital Conditions

% Known Initial Conditions for the ODE Solver
x1_0=0;
x2_0=0;
x3_0=0;
x4_0=0;
Lam_x1_0=0;
t0_0=0;

X0_Known_ODESolver=[x1_0;x2_0;x3_0;x4_0;Lam_x1_0];

% Initial Conditions from the Root Finder
Lam_x2_0=X0_RootFinder(1);
Lam_x3_0=X0_RootFinder(2);
Lam_x4_0=X0_RootFinder(3);
tf_0=X0_RootFinder(4);

X0_RootFinder_Guess=X0_RootFinder(1:end-1);

% Complete Initial Conditions for the ODE Solver
X0_ODESolver = [X0_Known_ODESolver;X0_RootFinder_Guess];

%% ODE Solver

% Creating Time Span for Integration
tf=tf_0;

T_Span_ODESolver=[0,tf];

% Options for the ODE Solver
Options=optimset('Display','Iter','TolFun',1e-6,'TolX',1e-6);

% Calling the ODE Solver
[t_ODESolver,X_sol_ODESolver]=ode113(@LinearTangent_ODEEquations ,...
    T_Span_ODESolver,X0_ODESolver,Options);

% Getting the Outputs from the ODE Solver
x1=X_sol_ODESolver(:,1);
x2=X_sol_ODESolver(:,2);
x3=X_sol_ODESolver(:,3);
x4=X_sol_ODESolver(:,4);
Lam_x1=X_sol_ODESolver(:,5);
Lam_x2=X_sol_ODESolver(:,6);
Lam_x3=X_sol_ODESolver(:,7);
Lam_x4=X_sol_ODESolver(:,8);

%% Computing the Hamiltonian for getting H_tf

% Getting Boundary Values
x1_tf=x1(end);
x2_tf=x2(end);
x3_tf=x3(end);
x4_tf=x4(end);
Lam_x1_tf=Lam_x1(end);
Lam_x2_tf=Lam_x2(end);
Lam_x3_tf=Lam_x3(end);
Lam_x4_tf=Lam_x4(end);

% Computing U - Control at tf
u_tf=atan2(Lam_x4_tf,Lam_x3_tf);

a=10; % Constant

% Computing H(tf)
H_tf=(Lam_x1_tf*x3_tf) + (Lam_x2_tf*x4_tf) + (Lam_x3_tf*(a*cos(u_tf))) +...
    (Lam_x4_tf*(a*sin(u_tf)));

%% Computing the Error

% Initializing Error
Error=zeros(5,1);

% Creating Error Vector
Error(1)=x2_tf-5;
Error(2)=x3_tf-45;
Error(3)=x4_tf-0;
Error(4)=Lam_x1_tf-0;
Error(5)=H_tf+1;

end

