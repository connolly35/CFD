%% CFD Boston University Course -------------------------------------------
% Step 7 of 12 to Navier Stokes Code: Two Dimmensional Diffusion
%
% Governing Equation: 
% \frac{\partial{u}}{\partial{t}} = \nu\left(\frac{\partial^2 u}{\partial x^2}+\frac{\partial^2 u}{\partial y^2}\right)
%
% This code is a second order spatial CFD problem, which is essentially the
% heat equation in two dimmensions. This is an explicit time marching - 
% upwind method with a FD in time (n) and a CD in space (i,j). The CD in 
% space is used because the diffusive property is isotropic
%
% Discretization is: 
% \frac{u^{n+1}_{i,j} - u^{n}_{i,j}}{\Delta t} = \frac{\nu}{\Delta x^2}(u^n_{i+1,j} -2u^n_{i,j} + u^n_{i-1,j}) + \frac{\nu}{\Delta y^2}(u^n_{i,j+1} -2u^n_{i,j} + u^n_{i,j-1})
%
% Solve Discretization for single unknown of the n+1 time variable and
% apply numberical method.
%
% Numerical Representation:
% u^{n+1}_{i,j} = u^n_{i,j} + \nu\frac{\Delta t}{\Delta x^2}(u^n_{i+1,j} - 2u^n_{i,j} + u^n_{i-1,j}) + \nu\frac{\Delta t}{\Delta y^2}(u^n_{i,j+1} - 2u^n_{i,j} + u^n_{i,j-1})
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
ny = nx;
x_final = 2;
y_final = 2;
dx = x_final/(nx - 1);
dy = y_final/(ny - 1);
x = 0:dx:x_final;
y = 0:dy:y_final;

%% Initial Conditions
for i = 1:nx
    for j = 1:ny
        if (x(i) >= 0.5 && x(i) <= 1) && (y(j) >= 0.5 && y(j) <= 1)
            u(i,j) = 2;
        else
            u(i,j) = 1;
        end
    end
end

%% Plot Initial Condition
figure(1);
plot(x,u,'--r');grid on;
xlabel('Axial Location')
ylabel('Velocity')
hold on;
% Three Dimensional Plot
figure(2);
h = mesh(x,y,u);
set(h,'LineStyle','--')
xlabel('Axial Location x')
ylabel('Axial Location y')
zlabel('Velocity')
hold on

%% Numerical Method (FD - Time, CD - Spacial)
for it = 1:nt
    un = u;
    for i = 2:nx-1
        for j = 2:nx-1
            u(i,j) = un(i,j) + nu*(dt/dx^2)*(un(i + 1,j) -2*un(i,j) + un(i - 1,j)) + nu*(dt/dy^2)*(un(i,j + 1) -2*un(i,j) + un(i,j - 1));
        end
    end
end

%% Plot Simulated Steady State Condition
figure(1); hold on;
plot(x,u,'b');grid on;
xlabel('Axial Location')
ylabel('Velocity')
% Three Dimensional Plot
figure(2);
mesh(x,y,u);
xlabel('Axial Location x')
ylabel('Axial Location y')
zlabel('Velocity')
hold on