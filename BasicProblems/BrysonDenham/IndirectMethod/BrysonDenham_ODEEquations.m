function [ Equation_Derivative ] = BrysonDenham_ODEEquations( T_Span_ODESolver,X0_ODESolver )

%% ODE Equations : Problem 1 - Part 2

%% Getting Required Values from the incoming Vectors

% From T_Span_ODESolver
t=T_Span_ODESolver;

% From X0_ODESolver
x=X0_ODESolver(1);
v=X0_ODESolver(2);
Lam_x=X0_ODESolver(3);
Lam_v=X0_ODESolver(4);

% Initializing P_Dot
Equation_Derivative=zeros(4,1);

%% Setting up the ODE Equations

% Computing U - Control
% u=-Lam_v;
u=-Lam_v;

% Equations
x_Derivative=v;
v_Derivative=u;
Lam_x_Derivative=0;
Lam_v_Derivative=-Lam_x;

% Creating Equation Vector
Equation_Derivative(1)=x_Derivative;
Equation_Derivative(2)=v_Derivative;
Equation_Derivative(3)=Lam_x_Derivative;
Equation_Derivative(4)=Lam_v_Derivative;

end