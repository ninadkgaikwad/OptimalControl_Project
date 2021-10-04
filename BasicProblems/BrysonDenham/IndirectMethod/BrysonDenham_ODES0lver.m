function [ Error ] = BrysonDenham_ODES0lver( X0_RootFinder )

%% ODE Solver and Error Calculator : Problem 1 - Part 2

%% Inital Conditions

% Known Initial Conditions for the ODE Solver
x_0=0;
v_0=1;
tf_0=1;

X0_Known_ODESolver=[x_0;v_0];

% Initial Conditions from the Root Finder
Lam_x_0=X0_RootFinder(1);
Lam_v_0=X0_RootFinder(2);

X0_RootFinder_Guess=X0_RootFinder;

% Complete Initial Conditions for the ODE Solver
X0_ODESolver = [X0_Known_ODESolver;X0_RootFinder_Guess];

%% ODE Solver

% Creating Time Span for Integration
tf=tf_0;

T_Span_ODESolver=[0,tf];

% Options for the ODE Solver
Options=optimset('Display','Iter','TolFun',1e-6,'TolX',1e-6);

% Calling the ODE Solver
[t_ODESolver,X_sol_ODESolver]=ode113(@BrysonDenham_ODEEquations ,...
    T_Span_ODESolver,X0_ODESolver,Options);

% Getting the Outputs from the ODE Solver
x=X_sol_ODESolver(:,1);
v=X_sol_ODESolver(:,2);
Lam_x=X_sol_ODESolver(:,3);
Lam_v=X_sol_ODESolver(:,4);

%% Computing the Hamiltonian for getting H_tf

% Getting Boundary Values
x_tf=x(end);
v_tf=v(end);
Lam_x_tf=Lam_x(end);
Lam_v_tf=Lam_v(end);

%% Computing the Error

% Initializing Error
Error=zeros(2,1);

% Creating Error Vector
Error(1)=x_tf-0;
Error(2)=v_tf+1;

end

