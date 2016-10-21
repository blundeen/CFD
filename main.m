clear; clc; close all; format compact;
%{
    AAE 5615
    Project 2
    Cody Webster
    Brandon Lundeen
%}


%% Parameters
gamma   = 1.4;      % [-]
R       = .287;     % [kJ/kg-K]
L       = 3;        % [m]
N       = 31;       % [-]
dx      = L/(N-1);  % [m]
n       = 1400;     % [iterations]

% Formula for areas
A       = @(x) 1+2.2*(x-1.5).^2;


%% Analytical Solution
% Declare Vectors for Non-Dimensional Variables
sol_M       = zeros(N,2);   % Possible Mach Numbers
anlyt.M     = zeros(N,1);   % Mach Numbers
anlyt.P     = zeros(N,1);   % Pressure [-]
anlyt.rho   = zeros(N,1);   % Density [-]
anlyt.T     = zeros(N,1);   % Temperature [-]
x           = 0:dx:L;       % Positions [m]
anlyt.A     = A(x);         % Analytical Areas [-]

syms M
for i = 1:length(x)
    eqn       	= (1/M^2)*((2/(gamma+1))*(1+((gamma-1)/2)*M^2))^((gamma+1)/(gamma-1))-anlyt.A(i)^2;
    sol     	  = solve(eqn,M);
    sol_M(i,1)	= abs(double(sol(1)));
    sol_M(i,2)	= abs(double(sol(2)));
    if x(i) <= 1.5
        anlyt.M(i)	= abs(double(sol(2)));
    else
        anlyt.M(i)	= abs(double(sol(1)));
    end
end
anlyt.P     = (1+((gamma-1)/2)*(anlyt.M).^2).^(-gamma/(gamma-1));
anlyt.rho   = (1+((gamma-1)/2)*(anlyt.M).^2).^(-1/(gamma-1));
anlyt.T     = (1+((gamma-1)/2)*(anlyt.M).^2).^(-1);

% Plot Analytical Data
figure('color',[1 1 1],'position',[680 250 560 730])
ax(1) = subplot(4,1,1);
plot(x,anlyt.M,'Linewidth',2)
ylabel('Mach Number','interpreter','latex')
grid on
ax(2) = subplot(4,1,2);
plot(x,anlyt.P,'Linewidth',2)
ylabel('$\displaystyle\frac{P}{P_0}$','interpreter','latex')
ylh = get(gca,'ylabel');    % Get Y-label Object Information
ylp = get(ylh, 'Position'); % Get Y-label Position Information
set(ylh, 'Rotation',0, 'Position',ylp, 'VerticalAlignment','middle', 'HorizontalAlignment','right')
grid on
ax(3) = subplot(4,1,3);
plot(x,anlyt.rho,'Linewidth',2)
ylabel('$\displaystyle\frac{\rho}{\rho_0}$','interpreter','latex')
ylh = get(gca,'ylabel');    % Get Y-label Object Information
ylp = get(ylh, 'Position'); % Get Y-label Position Information
set(ylh, 'Rotation',0, 'Position',ylp, 'VerticalAlignment','middle', 'HorizontalAlignment','right')
grid on
ax(4) = subplot(4,1,4);
plot(x,anlyt.T,'Linewidth',2)
xlabel('Position [m]','interpreter','latex')
ylabel('$\displaystyle\frac{T}{T_0}$','interpreter','latex')
ylh = get(gca,'ylabel');    % Get Y-label Object Information
ylp = get(ylh, 'Position'); % Get Y-label Position Information
set(ylh, 'Rotation',0, 'Position',ylp, 'VerticalAlignment','middle', 'HorizontalAlignment','right')
grid on
linkaxes(ax,'x')                                        % Link X-axis of all subplots
set(ax(1),'xticklabel',[],'position',[.13 .745 .8 .2])  % Remove X-tick values and setting plot positions
set(ax(2),'xticklabel',[],'position',[.13 .520 .8 .2])  % Remove X-tick values and setting plot positions
set(ax(3),'xticklabel',[],'position',[.13 .295 .8 .2])  % Remove X-tick values and setting plot positions
set(ax(4),'position',[.13 .06 .8 .2])                   % Setting plot positions

clearvars eqn sol ylh ylp ax

%% Experimental Solution




