%% CFD Boston University Course -------------------------------------------
% Step 4 of 12 to Navier Stokes Code: One Dimmensional Burgers Equation
%
% Governing Equation: 
% \frac{\partial{u}}{\partial{t}} + u\frac{\partial u}{\partial x} = \nu\frac{\partial^2 u}{\partial x^2}
%
% This code is a second order spatial CFD problem with convection, which is
% Burgers equation. This is an explicit time marching - upwind 
% method with a FD in time (n), BD in space (i) for the first derivative 
% and a CD in space (i) for the second derivative. The CD in space is
% used because the diffusive property is isotropic
%
% Discretization is: 
% \frac{u^{n+1}_{i} - u^{n}_{i}}{\Delta t} + u^n_i\frac{u^n_i - u^n_{i-1}}{\Delta x} = \frac{\nu}{\Delta x^2}(u^n_{i+1} -2u^n_i + u^n_{i-1})
%
% Solve Discretization for single unknown of the n+1 time variable and
% apply numberical method.
%
% Numerical Representation:
% u^{n+1}_i = u^n_i - u^n_i\frac{\Delta t}{\Delta x}(u^n_i - u^n_{i-1}) + \nu\frac{\Delta t}{\Delta x^2}(u^n_{i+1} -2u^n_i + u^n_{i-1})
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
dt = 0.000001;

%% Constants Setup
c = 1;
nu = 0.1;

%% Grid Setup
nx = 500;
x_final = 2;
dx = x_final/(nx - 1);
%x = 0:dx:x_final;
for i = 1:nx
    ip1(i) = i + 1;
    im1(i) = i - 1;
    x(i) = (i - 1)*dx;
end
ip1(nx) = 1;
im1(1) = nx;

%% Initial Conditions
for i = 1:nx
    phi = exp(-x(i)^2/4/nu) + exp(-(x(i)-2*pi).^2/4/nu);
    dphi = -0.5/nu*x(i)*exp(-x(i)^2/4/nu) - 0.5/nu*(x(i)-2*pi)*exp(-(x(i)-2*pi)^2/4/nu);
    u(i) = -2*nu*dphi/phi + 4;
end

%% Plot Initial Condition
figure(1);
plot(x,u,'--r');grid on;
xlabel('Axial Location')
ylabel('Velocity')
hold on;

%% Numerical Method (FD - Time, BD - Spacial First Derivative, and CD - Spacial Second Derivative)
% First Loop is in Time
for it = 1:nt
    t = (it - 1)*dt; % Time
    % Analytical Solution Loop
    for i = 1:nx
        phi = exp( -(x(i) -4*t)^2/4/nu/(t+1) ) + exp( -(x(i)-4*t-2*pi)^2/4/nu/(t+1) );
        dphi = -0.5/nu/(t+1)*( x(i)-4*t )*exp(-(x(i)-4*t)^2/4/nu/(t+1) ) - 0.5/nu/(t+1)*( x(i)-4*t-2*pi )*exp( -(x(i)-4*t-2*pi).^2/4/nu/(t+1) );
        ua(i) = -2*nu*dphi./phi + 4;
    end
    
    
    % Numerical Method Solution Loop
    un = u;
    for i = 2:nx-1
        u(i) = un(i) - un(i)*(dt/dx)*(un(i) - un(im1(i))) + nu*(dt/dx^2)*(un(ip1(i)) -2*un(i) + un(im1(i)));
    end
end

%% Plot Simulated Steady State Condition
figure(1); hold on;
plot(x,u,'b',x,ua,'-og');grid on;
xlabel('Axial Location')
ylabel('Velocity')
legend('Initial Condition','Numerical Solution','Analytical Solution');