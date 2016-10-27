clear; clc; close all; format compact;
%{
    AAE 5615
    Project 2
    Cody Webster
    Brandon Lundeen
%}

%% Parameters
gamma   = 1.4;      % [-]
R       = .287e3;   % [J/kg-K]
L       = 3;        % [m]
N       = 31;       % [-]
dx      = L/(N-1);  % [m]
n       = 1401;     % [iterations]

% Formula for areas
A       = @(x) 1+2.2*(x-1.5).^2;

%
%% Analytical Solution
% Declare Vectors for Non-Dimensional Variables
sol_M       = zeros(N,2);   % Possible Mach Numbers
anlyt.M     = zeros(N,1);   % Mach Numbers
anlyt.P     = zeros(N,1);   % Pressure [-]
anlyt.rho   = zeros(N,1);   % Density [-]
anlyt.T     = zeros(N,1);   % Temperature [-]
x           = 0:dx:L;       % Positions [m]
anlyt.x     = 0:dx:L;       % Positions [m]
anlyt.A     = A(x);         % Analytical Areas [-]

syms M
for i = 1:length(x)
    eqn       	= (1/M^2)*((2/(gamma+1))*(1+((gamma-1)/2)*M^2))^((gamma+1)/(gamma-1))-anlyt.A(i)^2;
    sol     	= solve(eqn,M);
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
%}
%
%% Iterative Solution for Part A-2
data_1.CFL      = .5;
x               = 0:dx:L;
q               = 1;

% Declares Vectors for data
data_1.rho      = zeros(n,N,q);
data_1.T        = zeros(n,N,q);
data_1.v        = zeros(n,N,q);

% Sets boundary conditions and initial conditions
data_1.rho(:,1,:)	= 1;
data_1.rho(1,:,:)   = 1-.3146*x;
data_1.T(:,1,:)     = 1;
data_1.T(1,:,:)     = 1-.2314*x;
data_1.v(1,:,:)     = (0.1+1.09*x).*(data_1.T(1,:,:).^.5);
data_1.A            = A(x);

% Declare vectors used in calculations
a_i         = zeros(N,1);
v_i         = zeros(N,1);
dt_i        = zeros(N-2,1);
data_1.dt   = zeros(n,1);
drho_n      = zeros(n,N,q);
dv_n        = zeros(n,N,q);
dT_n        = zeros(n,N,q);
rho_B       = zeros(n,N,q);
v_B         = zeros(n,N,q);
T_B         = zeros(n,N,q);
drho_B      = zeros(n,N,q);
dv_B        = zeros(n,N,q);
dT_B        = zeros(n,N,q);
drho_avg  	= zeros(n,N,q);
dv_avg     	= zeros(n,N,q);
dT_avg     	= zeros(n,N,q);

for i = 1:q             % 3rd Dimension of Matrix
    for j = 1:n         % 1st Dimension of Matrix
        % Calculates possible time steps to use for calculations
        for k = 2:N-1	% 2nd Dimension of Matrix
            a_i(k)          = sqrt(data_1.T(j,k,i));
            v_i(k)          = data_1.v(j,k,i);
            dt_i(j,k-1,i)	= data_1.CFL*(dx/(a_i(k)+v_i(k)));
        end
        
        % Determines actual time step to use for calculations
        data_1.dt(j)        = min(dt_i(j,:,i));
        
        % Predictor step of the calculations
        for k = 1:N-1	% 2nd Dimension of Matrix
            drho_n(j,k,i)       = -data_1.rho(j,k,i)*((data_1.v(j,k+1,i)-data_1.v(j,k,i))/dx)...
                -data_1.rho(j,k,i)*data_1.v(j,k,i)*((log(data_1.A(k+1))-log(data_1.A(k)))/dx)...
                -data_1.v(j,k,i)*((data_1.rho(j,k+1,i)-data_1.rho(j,k,i))/dx);
            dv_n(j,k,i)         = -data_1.v(j,k,i)*((data_1.v(j,k+1,i)-data_1.v(j,k,i))/dx)...
                -(1/gamma)*((data_1.T(j,k+1,i)-data_1.T(j,k,i))/dx...
                +(data_1.T(j,k,i)/data_1.rho(j,k,i))*(data_1.rho(j,k+1,i)-data_1.rho(j,k,i))/dx);
            dT_n(j,k,i)         = -data_1.v(j,k,i)*((data_1.T(j,k+1,i)-data_1.T(j,k,i))/dx)...
                -(gamma-1)*data_1.T(j,k,i)*((data_1.v(j,k+1,i)-data_1.v(j,k,i))/dx...
                +data_1.v(j,k,i)*((log(data_1.A(k+1))-log(data_1.A(k)))/dx));
            % ---------------------------------------
            rho_B(j,k,i)        = data_1.rho(j,k,i)+drho_n(j,k,i)*data_1.dt(j);
            v_B(j,k,i)          = data_1.v(j,k,i)  +dv_n(j,k,i)  *data_1.dt(j);
            T_B(j,k,i)          = data_1.T(j,k,i)  +dT_n(j,k,i)  *data_1.dt(j);
        end
        
        % Corrector Step of the Calculations
        for k = 2:N-1	% 2nd Dimension of Matrix
            drho_B(j,k,i)       = -rho_B(j,k,i)*((v_B(j,k,i)-v_B(j,k-1,i))/dx)...
                -rho_B(j,k,i)*v_B(j,k,i)*((log(data_1.A(k))-log(data_1.A(k-1)))/dx)...
                -v_B(j,k,i)*((rho_B(j,k,i)-rho_B(j,k-1,i))/dx);
            dv_B(j,k,i)         = -v_B(j,k,i)*((v_B(j,k,i)-v_B(j,k-1,i))/dx)...
                -(1/gamma)*(((T_B(j,k,i)-T_B(j,k-1,i))/dx)...
                +(T_B(j,k,i)/rho_B(j,k,i))*((rho_B(j,k,i)-rho_B(j,k-1,i))/dx));
            dT_B(j,k,i)         = -v_B(j,k,i)*((T_B(j,k,i)-T_B(j,k-1,i))/dx)...
                -(gamma-1)*T_B(j,k,i)*(((v_B(j,k,i)-v_B(j,k-1,i))/dx)...
                +v_B(j,k,i)*((log(data_1.A(k))-log(data_1.A(k-1)))/dx));
            % ---------------------------------------
            drho_avg(j,k,i)     = .5*(drho_n(j,k,i)+drho_B(j,k,i));
            dv_avg(j,k,i)       = .5*(dv_n(j,k,i)  +dv_B(j,k,i));
            dT_avg(j,k,i)       = .5*(dT_n(j,k,i)  +dT_B(j,k,i));
            % Data at Next Time Step ----------------
            data_1.rho(j+1,k,i) = data_1.rho(j,k,i)+drho_avg(j,k,i)*data_1.dt(j);
            data_1.v(j+1,k,i)	= data_1.v(j,k,i)  +dv_avg(j,k,i)  *data_1.dt(j);
            data_1.T(j+1,k,i)	= data_1.T(j,k,i)  +dT_avg(j,k,i)  *data_1.dt(j);
        end
        
        % Re-evaluate Boundary Conditions
        data_1.v(j+1,1,i)       = 2*data_1.v(j+1,2,i)    -data_1.v(j+1,3,i);
        data_1.rho(j+1,N,i)     = 2*data_1.rho(j+1,N-1,i)-data_1.rho(j+1,N-2,i);
        data_1.T(j+1,N,i)       = 2*data_1.T(j+1,N-1,i)  -data_1.T(j+1,N-2,i);
        data_1.v(j+1,N,i)       = 2*data_1.v(j+1,N-1,i)  -data_1.v(j+1,N-2,i);
    end
end

data_1.x = x;
data_1.P = data_1.rho.*data_1.T;            % Pressure [-]
data_1.M = data_1.v./sqrt(data_1.T);        % Mach Number [-]
data_1.dm = data_1.rho.*data_1.v.*data_1.A; % mass flow rate [-]

%% Part A-3
k = 16;
% Plot Analytical Data
figure('color',[1 1 1],'position',[680 250 560 730])
ax(1) = subplot(4,1,1);
plot(data_1.rho(:,k,1),'Linewidth',2)
ylabel('$\displaystyle\frac{\rho}{\rho_0}$','interpreter','latex')
ylh = get(gca,'ylabel');    % Get Y-label Object Information
ylp = get(ylh, 'Position'); % Get Y-label Position Information
set(ylh, 'Rotation',0, 'Position',ylp, 'VerticalAlignment','middle', 'HorizontalAlignment','right')
grid on
ax(2) = subplot(4,1,2);
plot(data_1.T(:,k,1),'Linewidth',2)
ylabel('$\displaystyle\frac{T}{T_0}$','interpreter','latex')
ylh = get(gca,'ylabel');    % Get Y-label Object Information
ylp = get(ylh, 'Position'); % Get Y-label Position Information
set(ylh, 'Rotation',0, 'Position',ylp, 'VerticalAlignment','middle', 'HorizontalAlignment','right')
grid on
ax(3) = subplot(4,1,3);
plot(data_1.P(:,k,1),'Linewidth',2)
ylabel('$\displaystyle\frac{P}{P_0}$','interpreter','latex')
ylh = get(gca,'ylabel');    % Get Y-label Object Information
ylp = get(ylh, 'Position'); % Get Y-label Position Information
set(ylh, 'Rotation',0, 'Position',ylp, 'VerticalAlignment','middle', 'HorizontalAlignment','right')
grid on
ax(4) = subplot(4,1,4);
plot(data_1.M(:,k,1),'Linewidth',2)
xlabel('Iteration','interpreter','latex')
ylabel('Mach Number','interpreter','latex')
grid on
linkaxes(ax,'x');                                       % Link X-axis of all subplots
set(ax(1),'xticklabel',[],'position',[.13 .745 .8 .2])  % Remove X-tick values and setting plot positions
set(ax(2),'xticklabel',[],'position',[.13 .520 .8 .2])  % Remove X-tick values and setting plot positions
set(ax(3),'xticklabel',[],'position',[.13 .295 .8 .2])  % Remove X-tick values and setting plot positions
set(ax(4),'position',[.13 .06 .8 .2])                   % Setting plot positions

%% Part A-4
figure
semilogy(abs(drho_avg(:,k,1)),'Linewidth',2)
xlabel('Iteration','interpreter','latex')
ylabel('$(\displaystyle\frac{d\rho}{dt})^{avg}$','interpreter','latex')
ylh = get(gca,'ylabel');    % Get Y-label Object Information
ylp = get(ylh, 'Position'); % Get Y-label Position Information
set(ylh, 'Rotation',0, 'VerticalAlignment','middle', 'HorizontalAlignment','center')
grid on

%% Part A-5
figure
ax(1) = subplot(2,1,1);
plot(x,data_1.rho(n,:,1),'Linewidth',2)
hold all
plot(x,anlyt.rho,'--','Linewidth',2)
hold off
ylabel('$\displaystyle\frac{\rho}{\rho_0}$','interpreter','latex')
ylh = get(gca,'ylabel');    % Get Y-label Object Information
ylp = get(ylh, 'Position'); % Get Y-label Position Information
set(ylh, 'Rotation',0, 'Position',ylp, 'VerticalAlignment','middle', 'HorizontalAlignment','right')
legend('Iterative','Analytical')
grid on
ax(2) = subplot(2,1,2);
plot(x,data_1.M(n,:,1),'Linewidth',2)
hold all
plot(x,anlyt.M,'--','Linewidth',2)
hold off
xlabel('Position [m]','interpreter','latex')
ylabel('Mach Number','interpreter','latex')
grid on
linkaxes(ax,'x');                                       % Link X-axis of all subplots
set(ax(1),'xticklabel',[],'position',[.13 .55 .8 .35])  % Remove X-tick values and setting plot positions
set(ax(2),'position',[.13 .15 .8 .35])                  % Setting plot positions

%% Part A-6
figure
plot(x,data_1.dm(1,:,1),'Linewidth',2)
hold all
plot(x,data_1.dm(101,:,1),'Linewidth',2)
plot(x,data_1.dm(701,:,1),'Linewidth',2)
plot(x,data_1.dm(1401,:,1),'Linewidth',2)
hold off
xlabel('Position [m]','interpreter','latex')
ylabel('Mass Flow Rate','interpreter','latex')
legend('i=1','i=100','i=700','i=1400')
grid on

%% Iterative Solution for Part A-7
data_2.CFL      = .5;
N               = 61;
dx              = L/(N-1);
x               = 0:dx:L;
q               = 1;

% Declares Vectors for data
data_2.rho      = zeros(n,N,q);
data_2.T        = zeros(n,N,q);
data_2.v        = zeros(n,N,q);

% Sets boundary conditions and initial conditions
data_2.rho(:,1,:)	= 1;
data_2.rho(1,:,:)   = 1-.3146*x;
data_2.T(:,1,:)     = 1;
data_2.T(1,:,:)     = 1-.2314*x;
data_2.v(1,:,:)     = (0.1+1.09*x).*(data_2.T(1,:,:).^.5);
data_2.A            = A(x);

% Declare vectors used in calculations
a_i         = zeros(N,1);
v_i         = zeros(N,1);
dt_i        = zeros(N-2,1);
data_2.dt   = zeros(n,1);
drho_n      = zeros(n,N,q);
dv_n        = zeros(n,N,q);
dT_n        = zeros(n,N,q);
rho_B       = zeros(n,N,q);
v_B         = zeros(n,N,q);
T_B         = zeros(n,N,q);
drho_B      = zeros(n,N,q);
dv_B        = zeros(n,N,q);
dT_B        = zeros(n,N,q);
drho_avg  	= zeros(n,N,q);
dv_avg     	= zeros(n,N,q);
dT_avg     	= zeros(n,N,q);

for i = 1:q             % 3rd Dimension of Matrix
    for j = 1:n         % 1st Dimension of Matrix
        % Calculates possible time steps to use for calculations
        for k = 2:N-1	% 2nd Dimension of Matrix
            a_i(k)          = sqrt(data_2.T(j,k,i));
            v_i(k)          = data_2.v(j,k,i);
            dt_i(j,k-1,i)	= data_2.CFL*(dx/(a_i(k)+v_i(k)));
        end
        
        % Determines actual time step to use for calculations
        data_2.dt(j)        = min(dt_i(j,:,i));
        
        % Predictor step of the calculations
        for k = 1:N-1	% 2nd Dimension of Matrix
            drho_n(j,k,i)       = -data_2.rho(j,k,i)*((data_2.v(j,k+1,i)-data_2.v(j,k,i))/dx)...
                -data_2.rho(j,k,i)*data_2.v(j,k,i)*((log(data_2.A(k+1))-log(data_2.A(k)))/dx)...
                -data_2.v(j,k,i)*((data_2.rho(j,k+1,i)-data_2.rho(j,k,i))/dx);
            dv_n(j,k,i)         = -data_2.v(j,k,i)*((data_2.v(j,k+1,i)-data_2.v(j,k,i))/dx)...
                -(1/gamma)*((data_2.T(j,k+1,i)-data_2.T(j,k,i))/dx...
                +(data_2.T(j,k,i)/data_2.rho(j,k,i))*(data_2.rho(j,k+1,i)-data_2.rho(j,k,i))/dx);
            dT_n(j,k,i)         = -data_2.v(j,k,i)*((data_2.T(j,k+1,i)-data_2.T(j,k,i))/dx)...
                -(gamma-1)*data_2.T(j,k,i)*((data_2.v(j,k+1,i)-data_2.v(j,k,i))/dx...
                +data_2.v(j,k,i)*((log(data_2.A(k+1))-log(data_2.A(k)))/dx));
            % ---------------------------------------
            rho_B(j,k,i)        = data_2.rho(j,k,i)+drho_n(j,k,i)*data_2.dt(j);
            v_B(j,k,i)          = data_2.v(j,k,i)  +dv_n(j,k,i)  *data_2.dt(j);
            T_B(j,k,i)          = data_2.T(j,k,i)  +dT_n(j,k,i)  *data_2.dt(j);
        end
        
        % Corrector Step of the Calculations
        for k = 2:N-1	% 2nd Dimension of Matrix
            drho_B(j,k,i)       = -rho_B(j,k,i)*((v_B(j,k,i)-v_B(j,k-1,i))/dx)...
                -rho_B(j,k,i)*v_B(j,k,i)*((log(data_2.A(k))-log(data_2.A(k-1)))/dx)...
                -v_B(j,k,i)*((rho_B(j,k,i)-rho_B(j,k-1,i))/dx);
            dv_B(j,k,i)         = -v_B(j,k,i)*((v_B(j,k,i)-v_B(j,k-1,i))/dx)...
                -(1/gamma)*(((T_B(j,k,i)-T_B(j,k-1,i))/dx)...
                +(T_B(j,k,i)/rho_B(j,k,i))*((rho_B(j,k,i)-rho_B(j,k-1,i))/dx));
            dT_B(j,k,i)         = -v_B(j,k,i)*((T_B(j,k,i)-T_B(j,k-1,i))/dx)...
                -(gamma-1)*T_B(j,k,i)*(((v_B(j,k,i)-v_B(j,k-1,i))/dx)...
                +v_B(j,k,i)*((log(data_2.A(k))-log(data_2.A(k-1)))/dx));
            % ---------------------------------------
            drho_avg(j,k,i)     = .5*(drho_n(j,k,i)+drho_B(j,k,i));
            dv_avg(j,k,i)       = .5*(dv_n(j,k,i)  +dv_B(j,k,i));
            dT_avg(j,k,i)       = .5*(dT_n(j,k,i)  +dT_B(j,k,i));
            % Data at Next Time Step ----------------
            data_2.rho(j+1,k,i) = data_2.rho(j,k,i)+drho_avg(j,k,i)*data_2.dt(j);
            data_2.v(j+1,k,i)	= data_2.v(j,k,i)  +dv_avg(j,k,i)  *data_2.dt(j);
            data_2.T(j+1,k,i)	= data_2.T(j,k,i)  +dT_avg(j,k,i)  *data_2.dt(j);
        end
        
        % Re-evaluate Boundary Conditions
        data_2.v(j+1,1,i)       = 2*data_2.v(j+1,2,i)    -data_2.v(j+1,3,i);
        data_2.rho(j+1,N,i)     = 2*data_2.rho(j+1,N-1,i)-data_2.rho(j+1,N-2,i);
        data_2.T(j+1,N,i)       = 2*data_2.T(j+1,N-1,i)  -data_2.T(j+1,N-2,i);
        data_2.v(j+1,N,i)       = 2*data_2.v(j+1,N-1,i)  -data_2.v(j+1,N-2,i);
    end
end

data_2.x = x;
data_2.P = data_2.rho.*data_2.T;            % Pressure [-]
data_2.M = data_2.v./sqrt(data_2.T);        % Mach Number [-]
data_2.dm = data_2.rho.*data_2.v.*data_2.A; % mass flow rate [-]

%% Part A-7
figure
ax(1) = subplot(2,1,1);
plot(data_1.x,data_1.rho(n,:,1),'Linewidth',2)
hold all
plot(data_2.x,data_2.rho(n,:,1),'Linewidth',2)
hold off
ylabel('$\displaystyle\frac{\rho}{\rho_0}$','interpreter','latex')
ylh = get(gca,'ylabel');    % Get Y-label Object Information
ylp = get(ylh, 'Position'); % Get Y-label Position Information
set(ylh, 'Rotation',0, 'Position',ylp, 'VerticalAlignment','middle', 'HorizontalAlignment','right')
legend('N=31','N=61')
grid on
ax(2) = subplot(2,1,2);
plot(data_1.x,data_1.M(n,:,1),'Linewidth',2)
hold all
plot(data_2.x,data_2.M(n,:,1),'Linewidth',2)
hold off
xlabel('Position [m]','interpreter','latex')
ylabel('Mach Number','interpreter','latex')
grid on
linkaxes(ax,'x');                                       % Link X-axis of all subplots
set(ax(1),'xticklabel',[],'position',[.13 .55 .8 .35])  % Remove X-tick values and setting plot positions
set(ax(2),'position',[.13 .15 .8 .35])                  % Setting plot positions
%}
%% Iterative Solution for Part A-8
data_3.CFL      = [1 1.5];
N               = 31;
dx              = L/(N-1);
x               = 0:dx:L;
q               = length(data_3.CFL);

% Declares Vectors for data
data_3.rho      = zeros(n,N,q);
data_3.T        = zeros(n,N,q);
data_3.v        = zeros(n,N,q);

% Sets boundary conditions and initial conditions
for i = 1:q
    data_3.rho(:,1,i)	= 1;
    data_3.rho(1,:,i)   = 1-.3146*x;
    data_3.T(:,1,i)     = 1;
    data_3.T(1,:,i)     = 1-.2314*x;
    data_3.v(1,:,i)     = (0.1+1.09*x).*(data_3.T(1,:,i).^.5);
    data_3.A            = A(x);
end

% Declare vectors used in calculations
a_i         = zeros(N,1);
v_i         = zeros(N,1);
dt_i        = zeros(N,1);
data_3.dt   = zeros(n,1,q);
drho_n      = zeros(n,N,q);
dv_n        = zeros(n,N,q);
dT_n        = zeros(n,N,q);
rho_B       = zeros(n,N,q);
v_B         = zeros(n,N,q);
T_B         = zeros(n,N,q);
drho_B      = zeros(n,N,q);
dv_B        = zeros(n,N,q);
dT_B        = zeros(n,N,q);
drho_avg  	= zeros(n,N,q);
dv_avg     	= zeros(n,N,q);
dT_avg     	= zeros(n,N,q);

for i = 1:q             % 3rd Dimension of Matrix
    for j = 1:n         % 1st Dimension of Matrix
        % Calculates possible time steps to use for calculations
%         for k = 2:N-1	% 2nd Dimension of Matrix
%             a_i(k)          = sqrt(data_3.T(j,k,i));
%             v_i(k)          = data_3.v(j,k,i);
%             dt_i(j,k-1,i)	= data_3.CFL(i)*(dx/(a_i(k)+v_i(k)));
%         end
        a_i(:,1)            = sqrt(data_3.T(j,:,i));
        v_i(:,1)            = data_3.v(j,:,i);
        dt_i(:,1)         = data_3.CFL(i).*(dx./(a_i(:,1)+v_i(:,1)));
        
        % Determines actual time step to use for calculations
        data_3.dt(j,1,i)        = min(dt_i(:,1));
        
        
        % Predictor step of the calculations
        for k = 1:N-1	% 2nd Dimension of Matrix
            drho_n(j,k,i)       = -data_3.rho(j,k,i)*((data_3.v(j,k+1,i)-data_3.v(j,k,i))/dx)...
                -data_3.rho(j,k,i)*data_3.v(j,k,i)*((log(data_3.A(k+1))-log(data_3.A(k)))/dx)...
                -data_3.v(j,k,i)*((data_3.rho(j,k+1,i)-data_3.rho(j,k,i))/dx);
            dv_n(j,k,i)         = -data_3.v(j,k,i)*((data_3.v(j,k+1,i)-data_3.v(j,k,i))/dx)...
                -(1/gamma)*((data_3.T(j,k+1,i)-data_3.T(j,k,i))/dx...
                +(data_3.T(j,k,i)/data_3.rho(j,k,i))*(data_3.rho(j,k+1,i)-data_3.rho(j,k,i))/dx);
            dT_n(j,k,i)         = -data_3.v(j,k,i)*((data_3.T(j,k+1,i)-data_3.T(j,k,i))/dx)...
                -(gamma-1)*data_3.T(j,k,i)*((data_3.v(j,k+1,i)-data_3.v(j,k,i))/dx...
                +data_3.v(j,k,i)*((log(data_3.A(k+1))-log(data_3.A(k)))/dx));
            % ---------------------------------------
            rho_B(j,k,i)        = data_3.rho(j,k,i)+drho_n(j,k,i)*data_3.dt(j,1,i);
            v_B(j,k,i)          = data_3.v(j,k,i)  +dv_n(j,k,i)  *data_3.dt(j,1,i);
            T_B(j,k,i)          = data_3.T(j,k,i)  +dT_n(j,k,i)  *data_3.dt(j,1,i);
        end
        
        % Corrector Step of the Calculations
        for k = 2:N-1	% 2nd Dimension of Matrix
            drho_B(j,k,i)       = -rho_B(j,k,i)*((v_B(j,k,i)-v_B(j,k-1,i))/dx)...
                -rho_B(j,k,i)*v_B(j,k,i)*((log(data_3.A(k))-log(data_3.A(k-1)))/dx)...
                -v_B(j,k,i)*((rho_B(j,k,i)-rho_B(j,k-1,i))/dx);
            dv_B(j,k,i)         = -v_B(j,k,i)*((v_B(j,k,i)-v_B(j,k-1,i))/dx)...
                -(1/gamma)*(((T_B(j,k,i)-T_B(j,k-1,i))/dx)...
                +(T_B(j,k,i)/rho_B(j,k,i))*((rho_B(j,k,i)-rho_B(j,k-1,i))/dx));
            dT_B(j,k,i)         = -v_B(j,k,i)*((T_B(j,k,i)-T_B(j,k-1,i))/dx)...
                -(gamma-1)*T_B(j,k,i)*(((v_B(j,k,i)-v_B(j,k-1,i))/dx)...
                +v_B(j,k,i)*((log(data_3.A(k))-log(data_3.A(k-1)))/dx));
            % ---------------------------------------
            drho_avg(j,k,i)     = .5*(drho_n(j,k,i)+drho_B(j,k,i));
            dv_avg(j,k,i)       = .5*(dv_n(j,k,i)  +dv_B(j,k,i));
            dT_avg(j,k,i)       = .5*(dT_n(j,k,i)  +dT_B(j,k,i));
            % Data at Next Time Step ----------------
            data_3.rho(j+1,k,i) = data_3.rho(j,k,i)+drho_avg(j,k,i)*data_3.dt(j,1,i);
            data_3.v(j+1,k,i)	= data_3.v(j,k,i)  +dv_avg(j,k,i)  *data_3.dt(j,1,i);
            data_3.T(j+1,k,i)	= data_3.T(j,k,i)  +dT_avg(j,k,i)  *data_3.dt(j,1,i);
        end
        
        % Re-evaluate Boundary Conditions
        data_3.v(j+1,1,i)       = 2*data_3.v(j+1,2,i)    -data_3.v(j+1,3,i);
        data_3.rho(j+1,N,i)     = 2*data_3.rho(j+1,N-1,i)-data_3.rho(j+1,N-2,i);
        data_3.T(j+1,N,i)       = 2*data_3.T(j+1,N-1,i)  -data_3.T(j+1,N-2,i);
        data_3.v(j+1,N,i)       = 2*data_3.v(j+1,N-1,i)  -data_3.v(j+1,N-2,i);
    end
end

data_3.x = x;
data_3.P = data_3.rho.*data_3.T;            % Pressure [-]
data_3.M = data_3.v./sqrt(data_3.T);        % Mach Number [-]
data_3.dm = data_3.rho.*data_3.v.*data_3.A; % mass flow rate [-]
close all
%% Part A-8
figure
ax(1) = subplot(2,1,1);
plot(data_1.x,data_1.rho(n,:,1),'Linewidth',2)
hold all
plot(data_3.x,data_3.rho(n,:,1),'Linewidth',2)
plot(data_3.x,data_3.rho(n,:,2),'Linewidth',2)
hold off
ylabel('$\displaystyle\frac{\rho}{\rho_0}$','interpreter','latex')
ylh = get(gca,'ylabel');    % Get Y-label Object Information
ylp = get(ylh, 'Position'); % Get Y-label Position Information
set(ylh, 'Rotation',0, 'Position',ylp, 'VerticalAlignment','middle', 'HorizontalAlignment','right')
legend('CFL = 0.5','CFL = 1','CFL = 1.5')
grid on
ax(2) = subplot(2,1,2);
plot(data_1.x,data_1.M(n,:,1),'Linewidth',2)
hold all
plot(data_3.x,data_3.M(n,:,1),'Linewidth',2)
plot(data_3.x,data_3.M(n,:,2),'Linewidth',2)
hold off
xlabel('Position [m]','interpreter','latex')
ylabel('Mach Number','interpreter','latex')
grid on
linkaxes(ax,'x');                                       % Link X-axis of all subplots
set(ax(1),'xticklabel',[],'position',[.13 .55 .8 .35])  % Remove X-tick values and setting plot positions
set(ax(2),'position',[.13 .15 .8 .35])                  % Setting plot positions
