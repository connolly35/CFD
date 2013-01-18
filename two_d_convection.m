%% CFD Boston University Course -------------------------------------------
% Step 6 of 12 to Navier Stokes Code: Two Dimmensional Convection
%
% Governing Equation: 
% \frac{\partial{u}}{\partial{t}} + u\frac{\partial u}{\partial x} + v\frac{\partial u}{\partial y} = 0 \\
% \frac{\partial{v}}{\partial{t}} + u\frac{\partial v}{\partial x} + v\frac{\partial v}{\partial y} = 0
%
% This code is a 2D nonlinear CFD problem, which is essentially the
% invisid Burgers equation. This is an explicit time marching - upwind 
% method with a FD in time (n) and a BD in space (i,j), It requires an
% equation for the x and y direction.
%
% Discretization is: 
% \frac{u^{n+1}_{i,j} - u^{n}_{i,j}}{\Delta t} + u^n_{i,j}\frac{u^n_{i,j} - u^n_{i-1,j}}{\Delta x} + v^n_{i,j}\frac{u^n_{i,j} - u^n_{i,j-1}}{\Delta y} = 0 \\
% \frac{v^{n+1}_{i,j} - v^{n}_{i,j}}{\Delta t} + u^n_{i,j}\frac{v^n_{i,j} - v^n_{i-1,j}}{\Delta x} + v^n_{i,j}\frac{v^n_{i,j} - v^n_{i,j-1}}{\Delta y} = 0
%
% Solve Discretization for single unknown of the n+1 time variable and
% apply numberical method.
%
% Numerical Representation:
% u^{n+1}_{i,j} = u^n_{i,j} - u^n_{i,j}\frac{\Delta t}{\Delta x}(u^n_{i,j} - u^n_{i-1,j}) - v^n_{i,j}\frac{\Delta t}{\Delta y}(u^n_{i,j} - u^n_{i,j-1})\\
% v^{n+1}_{i,j} = v^n_{i,j} - u^n_{i,j}\frac{\Delta t}{\Delta x}(v^n_{i,j} - v^n_{i-1,j}) - v^n_{i,j}\frac{\Delta t}{\Delta y}(v^n_{i,j} - v^n_{i,j-1})
%
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
            v(i,j) = 2;
        else
            u(i,j) = 1;
            v(i,j) = 1;
        end
    end
end

%% Plot Initial Condition
figure(1);
plot(x,u,'--r');grid on;
xlabel('Axial Location')
ylabel('Velocity')
hold on;
% Three Dimensional Plot - x
figure(2);
h = mesh(x,y,u);
set(h,'LineStyle','--')
xlabel('Axial Location x')
ylabel('Axial Location y')
zlabel('Velocity')
hold on
figure(3);
plot(x,v,'--r');grid on;
xlabel('Axial Location')
ylabel('Velocity')
hold on;
% Three Dimensional Plot - y
figure(4);
h = mesh(x,y,v);
set(h,'LineStyle','--')
xlabel('Axial Location x')
ylabel('Axial Location y')
zlabel('Velocity')
hold on

%% Numerical Method (FD - Time, BD - Spacial First Derivative for x and y)
% First Loop is in Time
for it = 1:nt
    un = u;
    vn = v;
    for i = 2:nx-1
        for j = 2:ny-1
            u(i,j) = un(i,j) - un(i,j)*(dt/dx)*(un(i,j) - un(i - 1,j)) - vn(i,j)*(dt/dy)*(un(i,j) - un(i,j - 1));
            v(i,j) = vn(i,j) - un(i,j)*(dt/dx)*(vn(i,j) - vn(i - 1,j)) - vn(i,j)*(dt/dy)*(vn(i,j) - vn(i,j - 1));
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
figure(3); hold on;
plot(x,v,'b');grid on;
xlabel('Axial Location')
ylabel('Velocity')
legend('Initial Condition','Numerical Solution','Analytical Solution');
% Three Dimensional Plot - y
figure(4);
mesh(x,y,v);
xlabel('Axial Location x')
ylabel('Axial Location y')
zlabel('Velocity')
hold on