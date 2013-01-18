%% CFD Boston University Course -------------------------------------------
% Step 2 of 12 to Navier Stokes Code: One Dimmensional Covection
%
% Governing Equation: 
% \frac{\partial{u}}{\partial{t}} + u\frac{\partial u}{\partial x} = 0
%
% This code is the simplest nonlinear CFD problem, which is essentially the
% invisid Burgers equation. This is an explicit time marching - upwind 
% method with a FD in time (n) and a BD in space (i)
%
% Discretization is: 
% \frac{u^{n+1}_{i} - u^{n}_{i}}{\Delta t} + u^n_i\frac{u^n_i - u^n_{i-1}}{\Delta x} = 0
%
% Solve Discretization for single unknown of the n+1 time variable and
% apply numberical method.
%
% Numerical Representation:
% u^{n+1}_i = u^n_i - u^n_i\frac{\Delta t}{\Delta x}(u^n_i - u^n_{i-1})
%
% Inputs: Boundary Conditions, Initial Conditions, and Grid
%
% Outputs: Final Plots of Results
%
% Created By:           Joseph Connolly
% Created On:           11/27/2011
%
% Modified On:
%
% -------------------------------------------------------------------------

%% Clear Old Data
clear all;close all;

%% Simulation Time Setup
nt = 50;
dt = 0.001;

%% Constants Setup
c = 1;

%% Grid Setup
nx = 100;
x_final = 2;
dx = x_final/(nx - 1);
x = 0:dx:x_final;

%% Initial Conditions
for i = 1:nx
    if x(i) >= 0.5 && x(i) <= 1;
        u(i) = 2;
    else
        u(i) = 1;
    end
end

%% Plot Initial Condition
figure(1);
plot(x,u,'--r');grid on;
xlabel('Axial Location')
ylabel('Velocity')
hold on;

%% Numerical Method (FD - Time, BD - Spacial)
for it = 1:nt
    un = u;
    for i = 2:nx-1
        u(i) = un(i) - un(i)*(dt/dx)*(un(i) - un(i - 1));
    end
end

%% Plot Simulated Steady State Condition
figure(1); hold on;
plot(x,u,'b');grid on;
xlabel('Axial Location')
ylabel('Velocity')