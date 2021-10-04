function [ Equation_Derivative ] = LinearTangent_ODEEquations( T_Span_ODESolver,X0_ODESolver )

%% ODE Equations : Problem 1 - Part 2

%% Getting Required Values from the incoming Vectors

% From T_Span_ODESolver
t=T_Span_ODESolver;

% From X0_ODESolver
x1=X0_ODESolver(1);
x2=X0_ODESolver(2);
x3=X0_ODESolver(3);
x4=X0_ODESolver(4);
Lam_x1=X0_ODESolver(5);
Lam_x2=X0_ODESolver(6);
Lam_x3=X0_ODESolver(7);
Lam_x4=X0_ODESolver(8);

% Initializing P_Dot
Equation_Derivative=zeros(8,1);

%% Setting up the ODE Equations

% Constant
a=100;

% Computing U - Control
u=atan2(Lam_x4,Lam_x3);

% Equations
x1_Derivative=x3;
x2_Derivative=x4;
x3_Derivative=a*cos(u);
x4_Derivative=a*sin(u);
Lam_x1_Derivative=0;
Lam_x2_Derivative=0;
Lam_x3_Derivative=-Lam_x1;
Lam_x4_Derivative=-Lam_x2;

% Creating Equation Vector
Equation_Derivative(1)=x1_Derivative;
Equation_Derivative(2)=x2_Derivative;
Equation_Derivative(3)=x3_Derivative;
Equation_Derivative(4)=x4_Derivative;
Equation_Derivative(5)=Lam_x1_Derivative;
Equation_Derivative(6)=Lam_x2_Derivative;
Equation_Derivative(7)=Lam_x3_Derivative;
Equation_Derivative(8)=Lam_x4_Derivative;

end