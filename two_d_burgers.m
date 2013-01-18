%% CFD Boston University Course -------------------------------------------
% Step 8 of 12 to Navier Stokes Code: Two Dimmensional Burgers Equation
%
% Governing Equation: 
% \frac{\partial{u}}{\partial{t}} + u\frac{\partial u}{\partial x} + v\frac{\partial u}{\partial y} = \nu\left(\frac{\partial^2 u}{\partial x^2}+\frac{\partial^2 u}{\partial y^2}\right)
% \frac{\partial{v}}{\partial{t}} + u\frac{\partial v}{\partial x} + v\frac{\partial v}{\partial y} = \nu\left(\frac{\partial^2 v}{\partial x^2}+\frac{\partial^2 v}{\partial y^2}\right)
%
% This code is a second order spatial CFD problem with convection, which is
% Burgers 2D equation. This is an explicit time marching - upwind 
% method with a FD in time (n), BD in space (i,j) for the first derivative 
% and a CD in space (i,j) for the second derivative. The CD in space is
% used because the diffusive property is isotropic
%
% Discretization is: 
% \frac{u^{n+1}_{i,j} - u^{n}_{i,j}}{\Delta t} + u^n_{i,j}\frac{u^n_{i,j} - u^n_{i-1,j}}{\Delta x} + v^n_{i,j}\frac{u^n_{i,j} - u^n_{i,j-1}}{\Delta y} = {\nu}\left(\frac{u^n_{i+1,j} -2u^n_{i,j} + u^n_{i-1,j}}{\Delta x^2} + \frac{u^n_{i,j+1} -2u^n_{i,j} + u^n_{i,j-1}}{\Delta y^2}\right)\\
% \frac{v^{n+1}_{i,j} - v^{n}_{i,j}}{\Delta t} + u^n_{i,j}\frac{v^n_{i,j} - v^n_{i-1,j}}{\Delta x} + v^n_{i,j}\frac{v^n_{i,j} - v^n_{i,j-1}}{\Delta y} = {\nu}\left(\frac{v^n_{i+1,j} -2v^n_{i,j} + v^n_{i-1,j}}{\Delta x^2} + \frac{v^n_{i,j+1} -2v^n_{i,j} + v^n_{i,j-1}}{\Delta y^2}\right)
%
% Solve Discretization for single unknown of the n+1 time variable and
% apply numberical method for both x and y direction.
%
% Numerical Representation:
% u^{n+1}_{i,j} = u^n_{i,j} - u^n_{i,j}\frac{\Delta t}{\Delta x}(u^n_{i,j} - u^n_{i-1,j}) - v^n_{i,j}\frac{\Delta t}{\Delta y}(u^n_{i,j} - u^n_{i,j-1}) +\dots\\
% \dots\nu{\Delta t}\left(\frac{u^n_{i+1,j} -2u^n_{i,j} + u^n_{i-1,j}}{\Delta x^2}+\frac{u^n_{i,j+1} -2u^n_{i,j} + u^n_{i,j-1}}{\Delta y^2}\right) \\
% v^{n+1}_{i,j} = v^n_{i,j} - u^n_{i,j}\frac{\Delta t}{\Delta x}(v^n_{i,j} - v^n_{i-1,j}) - v^n_{i,j}\frac{\Delta t}{\Delta y}(v^n_{i,j} - v^n_{i,j-1}) +\dots\\
% \dots\nu{\Delta t}\left(\frac{v^n_{i+1,j} -2v^n_{i,j} + v^n_{i-1,j}}{\Delta x^2}+\frac{v^n_{i,j+1} -2v^n_{i,j} + v^n_{i,j-1}}{\Delta y^2}\right)
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

%% Numerical Method (FD - Time, BD - Spacial First Derivative, and CD - Spacial Second Derivative)
% First Loop is in Time
for it = 1:nt
    un = u;
    vn = v;
    for i = 2:nx-1
        for j = 2:ny-1
            u(i,j) = un(i,j) - un(i,j)*(dt/dx)*(un(i,j) - un(i - 1,j)) - vn(i,j)*(dt/dy)*(un(i,j) - un(i,j - 1)) + nu*(dt/dx^2)*(un(i + 1,j) -2*un(i,j) + un(i - 1,j)) + nu*(dt/dy^2)*(un(i,j + 1) -2*un(i,j) + un(i,j - 1));
            v(i,j) = vn(i,j) - un(i,j)*(dt/dx)*(vn(i,j) - vn(i - 1,j)) - vn(i,j)*(dt/dy)*(vn(i,j) - vn(i,j - 1)) + nu*(dt/dx^2)*(vn(i + 1,j) -2*vn(i,j) + vn(i - 1,j)) + nu*(dt/dy^2)*(vn(i,j + 1) -2*vn(i,j) + vn(i,j - 1));
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