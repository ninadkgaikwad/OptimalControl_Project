%% Main File : Problem 1- Part2 : Indirect Shooting

clear all;
clc;
close all;

%% Initial Conditions

% Known Initial Conditions for the ODE Solver
x_0=0;
v_0=1;
tf_0=1;

X0_Known_ODESolver=[x_0;v_0;tf_0];

% Unknown Initial Conditions for the Root Finder
% Lam_x_0=randn(1);
% Lam_v_0=randn(1);

Lam_x_0=randn(1);
Lam_v_0=[-6];

X0_RootFinder=[Lam_x_0;Lam_v_0];

%% Root Finder

% Options for the Root Finder
Options=optimset('Display','iter','TolFun',1e-6);

% Calling the Root Finder
[X_Sol_RootFinder,fval,exitflag,output] = ...
    fsolve(@BrysonDenham_ODES0lver,X0_RootFinder,Options);

%% Solving the ODE with Optimal Initial Conditions

% Inital Conditions

% Known Initial Conditions for the ODE Solver
x_0=0;
v_0=1;

% Initial Conditions from the Root Finder
Lam_x_0=X_Sol_RootFinder(1);
Lam_v_0=X_Sol_RootFinder(2);

% Lam_x_0=[-6300.95257935419];
% Lam_v_0=[-2046.13130763769];

% Complete Initial Conditions for the ODE Solver
X0_ODESolver = [x_0;v_0;Lam_x_0;Lam_v_0];

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

% Computing Control

% Control=-Lam_v;
Control=Lam_v;

%% Plotting Results

% Plotting States 
figure(1)
hold on
plot(t_ODESolver,x,'-g','LineWidth',1.5);
plot(t_ODESolver,v,'-.m','LineWidth',1.5);
title('States vs. Time','Interpreter','latex');
xlabel('Time (sec)','Interpreter','latex');
ylabel('States','Interpreter','latex');
legend1=legend('$x$','$v$');
set(legend1,'Interpreter','latex');
hold off;

% Plotting Control
figure(2)
hold on
plot(t_ODESolver,Control,'-k','LineWidth',1.5);
title('Control - $u(t)$ vs. Time','Interpreter','latex');
xlabel('Time (sec)','Interpreter','latex');
ylabel('$u(t)$','Interpreter','latex');
hold off;

% Plotting Costates
figure(3)
hold on
plot(t_ODESolver,Lam_x,'-g','LineWidth',1.5);
plot(t_ODESolver,Lam_v,'-.m','LineWidth',1.5);
title('Co-States vs. Time','Interpreter','latex');
xlabel('Time (sec)','Interpreter','latex');
ylabel('Co-States','Interpreter','latex');
legend2=legend('$\lambda_x$','$\lambda_v$');
set(legend2,'Interpreter','latex');
hold off;