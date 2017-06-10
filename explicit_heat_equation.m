% 1D Heat Flow: Flux and Conservation
% Harlow (1995)
% Bryan Kaiser
% 3/18/14

close all 
clear all
clc

% Definitions:
% Flux:the amount of something passing through a unit area in a unit time.
% (from p.7)

% Computational Set Up and Notes:
% - dx = D/j = zone length.
% - zones go from 0 to j+1 (halo zones included).
% - Temperatures defined at each cell center.
% - V = A*dx.
% - Heat energy of cell = m*b*T, where b = specific heat
% - m = rho*A*dx
% - Heat energy of cell = (rho*A*dx)*b*T
% - t(n) = n*dt
% - k/(b*rho) = thermometric conductivity of a material = sigma

% Governing Energy Conservation Equation:
% 1) (Heat E_j)^(n+1) - (Heat E_j)^(n) =
%                     = dt*A*[(Flux_j-1/2)^(n)-(Flux_j+1/2)^(n)]
% or: dE/dt = k*A*(dT/dx|_j-1/2 - dT/dx|_j+1/2)
%

% Definitions
% T_L = left end temperature, on the left edge of cell T2 (T1 in the book)
% T_L is constant, so T_L = (T1+T2)/2
% T_R is also constant, defined the same way
% T0 = 2*T_L-T1
% Tj+1 = 2*T_R-Tj

%==========================================================================

% Variables
T0 = 0; % C, Initial rod temperature
TL = 400; % C, Fixed left temperature
TR = 0; % C, Fixed right temperature
sigma = 1; % m^2, thermometric conductivity
dx = 1; % m, zone length
Nx = 50; % number of zones
dt = 0.1; % s, time step
time = 500; % s, total time
X = zeros(1,Nx);
X(1) = -dx;
for i = 2:Nx
    X(i) = X(i-1)+dx;
end

% Diffusion stability limit (for this problem):
DSL = sigma*dt/dx^2 % Should be <0.5

% Sections for FDM
% 1) Initialization
% 2) Loop set up for time cycle
% 3) Define BCs and update ghost zone values
% 4) Implement conservation equation (two loops, assign Tnew, transfer
%    values of Tnew back into T.
% 5) Output subroutine

%--------------------------------------------------------------------------

% 1) Initialization
T = zeros(2,Nx);
T(1,:) = T(1,:)+T0; % Initial condition at t=0
t = 0;
old = 1;
new = 2;
countm = 1;

% 2) Loop Set Up
while t <= (time-dt)
    % Time advancement
    t = t + dt;
    
    % 3) Define BCs and update ghost zone values
    T(old,1) = 2*TL-T(old,2); % BC left
    T(old,end) = 2*TR-T(old,end-1); % BC right
   
    % Movie Output
    figure(1)
    set(gca,'FontSize',14)
    plot(X,T(new,:),'r')
    axis([0,50,0,400]);
    grid on
    xlabel('x (m)')
    ylabel('T (C)')
    %title('Continuity')
    M(countm)=getframe;
    countm = countm+1;
     
    % 4) Implement conservation equation 
    % Assign T(new):
    for j = 2:(Nx-1)
     T(new,j) = T(old,j)+sigma*dt/(dx^2)*[T(old,j+1)-...
         2*T(old,j)+T(old,j-1)];
    end
    
    % Transfer values of T(new) back to T(old)
    T(old,:) = T(new,:);
    
end




