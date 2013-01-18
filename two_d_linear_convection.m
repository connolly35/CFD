%% CFD Boston University Course -------------------------------------------
% Step 5 of 12 to Navier Stokes Code: Two Dimmensional Linear Convection
%
% Governing Equation: 
% \frac{\partial{u}}{\partial{t}} + c\frac{\partial u}{\partial x} + c\frac{\partial u}{\partial y} = 0
%
% This code is the first simplest 2D CFD problem, which is essentially a 
% wave propagation problem. This is an explicit time marching - upwind 
% method with a FD in time (n) and a BD in space (i,j)
%
% Discretization is: 
% \frac{u^{n+1}_{i,j} - u^{n}_{i,j}}{\Delta t} + C\frac{u^n_{i,j} - u^n_{i-1,j}}{\Delta x} + C\frac{u^n_{i,j} - u^n_{i,j-1}}{\Delta y} = 0
%
% Solve Discretization for single unknown of the n+1 time variable and
% apply numberical method.
%
% Numerical Representation:
% u^{n+1}_{i,j} = u^n_{i,j} - c\frac{\Delta t}{\Delta x}(u^n_{i,j} - u^n_{i-1,j}) - c\frac{\Delta t}{\Delta y}(u^n_{i,j} - u^n_{i,j-1})
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
nt = 500;
dt = 0.0001;

%% Constants Setup
c = 1;
nu = 0.1;

%% Grid Setup
nx = 100;
ny = nx; % For Mesh Plots they Must be the Same Length
x_final = 2;
y_final = 2;
dx = x_final/(nx - 1);
dy = y_final/(ny - 1);
% Define 2D Mesh
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

%% Numerical Method (FD - Time, BD - Spacial First Derivative for x and y)
% First Loop is in Time
for it = 1:nt
    un = u;
    for i = 2:nx-1
        for j = 2:ny-1
            u(i,j) = un(i,j) - c*(dt/dx)*(un(i,j) - un(i - 1,j)) - c*(dt/dy)*(un(i,j) - un(i,j - 1));
        end
    end
end

%% Plot Simulated Steady State Condition
figure(1); hold on;
plot(x,u,'b');grid on;
xlabel('Axial Location')
ylabel('Velocity')
legend('Initial Condition','Numerical Solution','Analytical Solution');
% Three Dimensional Plot
figure(2);
mesh(x,y,u);
xlabel('Axial Location x')
ylabel('Axial Location y')
zlabel('Velocity')
hold on