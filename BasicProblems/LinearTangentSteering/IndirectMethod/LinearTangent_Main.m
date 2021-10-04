%% Main File : Linear Tangent Steering : Indirect Shooting

clear all;
clc;
close all;

%% Initial Conditions

% Known Initial Conditions for the ODE Solver
x1_0=0;
x2_0=0;
x3_0=0;
x4_0=0;
Lam_x1_0=0;
t0_0=0;

X0_Known_ODESolver=[x1_0;x2_0;x3_0;x4_0;Lam_x1_0;t0_0];

% Unknown Initial Conditions for the Root Finder
Lam_x2_0=randn(1);
Lam_x3_0=randn(1);
Lam_x4_0=randn(1);
tf_0=1;

X0_RootFinder=[Lam_x2_0;Lam_x3_0;Lam_x4_0;tf_0];

%% Root Finder

% Options for the Root Finder
Options=optimset('Display','iter','TolFun',1e-5);

% Calling the Root Finder
tic; % Timing the Process

[X_Sol_RootFinder,fval,exitflag,output] = ...
    fsolve(@LinearTangent_ODES0lver,X0_RootFinder,Options);

TimeTaken = toc;

%% Solving the ODE with Optimal Initial Conditions

% Inital Conditions

% Known Initial Conditions for the ODE Solver
x1_0=0;
x2_0=0;
x3_0=0;
x4_0=0;
Lam_x1_0=0;
t0_0=0;

% Initial Conditions from the Root Finder
Lam_x2_0=X_Sol_RootFinder(1);
Lam_x3_0=X_Sol_RootFinder(2);
Lam_x4_0=X_Sol_RootFinder(3);
tf_0=X_Sol_RootFinder(4);

% Complete Initial Conditions for the ODE Solver
X0_ODESolver = [x1_0;x2_0;x3_0;x4_0;Lam_x1_0;Lam_x2_0;Lam_x3_0;Lam_x4_0];

% Creating Time Span for Integration
tf=tf_0;
T_Span_ODESolver=[0,tf];

% Options for the ODE Solver
Options=optimset('Display','Iter','TolFun',1e-3,'TolX',1e-3);

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

% Computing Control

Control=atan2(Lam_x4,Lam_x3);

%% Plotting Results

% Plotting States 
figure(1)
hold on
grid on
plot(t_ODESolver,x1,'-g','LineWidth',1.5);
plot(t_ODESolver,x2,'-.m','LineWidth',1.5);
plot(t_ODESolver,x3,'-r','LineWidth',1.5);
plot(t_ODESolver,x4,'--b','LineWidth',1.5);
title('States vs. Time','Interpreter','latex');
xlabel('Time ','Interpreter','latex');
ylabel('States','Interpreter','latex');
legend1=legend('$x_{1}$','$x_{2}$','$x_{3}$','$x_{4}$');
set(legend1,'Interpreter','latex');
hold off;

% Plotting Control
figure(2)
hold on
grid on
plot(t_ODESolver,Control,'-k','LineWidth',1.5);
title('Control - $u(t)$ vs. Time','Interpreter','latex');
xlabel('Time ','Interpreter','latex');
ylabel('$u(t)$','Interpreter','latex');
hold off;

% Plotting Costates
figure(3)
hold on
grid on
plot(t_ODESolver,Lam_x1,'-g','LineWidth',1.5);
plot(t_ODESolver,Lam_x2,'-.m','LineWidth',1.5);
plot(t_ODESolver,Lam_x3,'-r','LineWidth',1.5);
plot(t_ODESolver,Lam_x4,'--b','LineWidth',1.5);
title('Co-States vs. Time','Interpreter','latex');
xlabel('Time ','Interpreter','latex');
ylabel('Co-States','Interpreter','latex');
legend2=legend('$\lambda_{x_{1}}$','$\lambda_{x_{2}}$','$\lambda_{x_{3}}$','$\lambda_{x_{4}}$');
set(legend2,'Interpreter','latex');
hold off;

fprintf('Time taken to solve the HBVP = %.4f',TimeTaken)