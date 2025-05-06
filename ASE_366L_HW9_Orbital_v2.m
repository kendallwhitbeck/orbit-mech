%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ASE 366L: Applied Orbital Mechanics
% HOMEWORK 9 MAIN FILE
tic
close all; format long g; clc; clear all;
% set(0, 'DefaultAxesFontSize', 10, 'defaultlinelinewidth',1)

%% Problem 2, Part (a.)
% Initial Relative Position & Velocity (RSW frame)
rho_0 = [10 0 0]'; % m
zeta_0 = [0 -.1 -.2]'; % m/s

sma = 6700; % Semi-Major Axis (km)
mu_earth = 398600.4415; % km^3/s^2

% Calc Mean Motion
n = sqrt(mu_earth/sma^3); % Hertz

% Final Time in Seconds
t_f = 45*60;

% Calc Final Relative Position in RSW frame
[rho_f,zeta_f] = CW_soln(rho_0,zeta_0,t_f,n);

fprintf('-------------- PROBLEM 2, PART (a) -------------- \n\n')
fprintf('The Final Relative Position (RSW frame) of \nthe chaser from the target in meters is:\n')
disp(rho_f)

%% Problem 2, Part (b.)
zeta_S_C = -6*n*rho_0(1) - 3 * zeta_0(2); % m/s
rho_R_C = -2*zeta_S_C / (3*n); % m
rho_S_C = rho_0(2) - 2*zeta_0(1)/n; % m
rho_W_C = rho_0(3); % m

C = sqrt((3*rho_0(1) + 2*zeta_0(2)/n)^2 + (zeta_0(1)/n)^2);

D = sqrt((zeta_0(3)/n)^2 + rho_0(3)^2);

fprintf('-------------- PROBLEM 2, PART (b) -------------- \n\n')
fprintf('rho_R_C in meters is:\n')
disp(rho_R_C)
fprintf('rho_S_C in meters is:\n')
disp(rho_S_C)
fprintf('rho_W_C in meters is:\n')
disp(rho_W_C)
fprintf('zeta_S_C (their change over time) in m/s is:\n')
disp(zeta_S_C)
fprintf('Axis Length C in meters is:\n')
disp(C)
fprintf('Axis Length D in meters is:\n')
disp(D)

%% Problem 5, Part (b.)
% Initial Relative: Position (rho_05), Velocity (zeta_05)
rho_05 = [25; 0; 0]; % m
zeta_05 = [-1; 0; 0]; % m/s

% Calculating necessary axes
R_unit = [rho_05 / norm(rho_05)];
W_unit = [cross(rho_05,zeta_05) ./ (norm(rho_05) * norm(zeta_05))];
S_unit = cross(W_unit,R_unit);

% Defining Transformation Matrix for RSW -> ijk
T_RSW_ijk = [R_unit S_unit W_unit];

% Calc chaser position & velocity in inertial frame
r = T_RSW_ijk * rho_05 ./ 1000 % km
v = T_RSW_ijk * zeta_05 ./ 1000 % km/s


%% Problem 5, Part (c.)
sma = 7000; % km
ecc = 0;
inc = 75*pi/180; % rad
RAAN = 0; % rad
argp = 0; % rad
tru = 0; % rad

% Vector of Keplarian Orbital Elements
xKOE = [sma, ecc, inc, RAAN, argp, tru];

% Convert initial KOEs into Cartesian position/velocity
Cart0 = kep2cart_Whitbeck(mu_earth,xKOE);

% % Calculate Orbit Period, T, & create matching timespan
% T = ceil(2 * pi * sqrt(sma^3 / mu_earth)); % sec
tspan = (0:1:t_f); % 1 orbit period (sec)
% [m,n] = size(tspan);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ODE45 Setup %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Anonymous function outputting 6x1 vector w/ 3D velocity & acceleration
% w/ input of time & 6x1 3D position/velocity vector.
vel_acc = @(t,Cart) [Cart(4:6);...
    -mu_earth .* Cart(1:3) ./ (norm(Cart(1:3))^3)];

% ODE options as given in lecture
odeoptions = odeset('RelTol',1e-28,'AbsTol',1e-30);

% Calling ode45 function to output time vector T & position/velocity matrix
% at each time interval Y.
[~,Y] = ode45(vel_acc, tspan, Cart0, odeoptions);

% Get Y matrix in form where columns contain pos/vel & each new column is
% another time t, then break it into pos matrix & vel matrix
Y = Y';









%%
time_to_run_in_seconds = toc
