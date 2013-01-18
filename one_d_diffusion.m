%% CFD Boston University Course -------------------------------------------
% Step 3 of 12 to Navier Stokes Code: One Dimmensional Diffusion
%
% Governing Equation: 
% \frac{\partial{u}}{\partial{t}} = \nu\frac{\partial^2 u}{\partial x^2}
%
% This code is a second order spatial CFD problem, which is essentially the
% heat equation. This is an explicit time marching - upwind 
% method with a FD in time (n) and a CD in space (i). The CD in space is
% used because the diffusive property is isotropic
%
% Discretization is: 
% \frac{u^{n+1}_{i} - u^{n}_{i}}{\Delta t} = \frac{\nu}{\Delta x^2}(u^n_{i+1} -2u^n_i + u^n_{i-1})
%
% Solve Discretization for single unknown of the n+1 time variable and
% apply numberical method.
%
% Numerical Representation:
% u^{n+1}_i = u^n_i + \nu\frac{\Delta t}{\Delta x^2}(u^n_{i+1} -2u^n_i + u^n_{i-1})
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
nu = 0.1;

%% Grid Setup
nx = 50;
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

%% Numerical Method (FD - Time, CD - Spacial)
for it = 1:nt
    un = u;
    for i = 2:nx-1
        u(i) = un(i) + nu*(dt/dx^2)*(un(i + 1) -2*un(i) + un(i - 1));
    end
end

%% Plot Simulated Steady State Condition
figure(1); hold on;
plot(x,u,'b');grid on;
xlabel('Axial Location')
ylabel('Velocity')