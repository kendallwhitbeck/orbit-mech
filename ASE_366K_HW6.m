clc; clear all; close all; format shortg;
set(0, 'DefaultAxesFontSize', 16, 'defaultlinelinewidth',1)

%% Problem 1
% Givens
mu = 398600.4415; % km^3/s^2
R = 6378.1363; % km
J2 = 0.00108248;

% Inertial, geocentric, ijk frame at time t0 = 0
% Initial Position vector
r0 = [-140.879 6715.125 1291.458]; % km
% Initial Velocity vector
v0 = [7.472003 0.511221 -1.4448890]; % km/s
% Initial Position/Velocity vector
X0 = [r0, v0]';
% Initial Energy scalar using Vis-Viva
energy0 = dot(v0,v0) / 2 - mu / norm(r0);


%% ODE45 Setup
% Anonymous function outputting 6x1 vector w/ 3D velocity & acceleration
% w/ input of time & 6x1 3D position/velocity vector.
vel_acc = @(t,r) [r(4:6); -mu .* r(1:3) ./ (norm(r(1:3))^3)];
    
% Time span every 20 seconds starting w/ initial time t0 & ending w/ final
% time tf.
t0 = 0; % seconds
tf = 16880; % seconds
tspan = (t0:20:tf);
tspanhr = tspan / 3600; % hours

% ODE options as given in lecture
odeoptions = odeset('RelTol',1e-13,'AbsTol',1e-20);

% Calling ode45 function to output time vector T & position/velocity matrix
% at each time interval Y.
[~,Y] = ode45(vel_acc, tspan, X0, odeoptions);

% Get Y matrix in form where columns contain pos/vel & each new column is
% another time t
Y = Y';

%% #1 Part (a)
% Calc magnitudes of position, velocity & acceleration vecs over time
n = length(tspan);
r_mag_vec = zeros(1,n);
v_mag_vec = zeros(1,n);
a_mag_vec = zeros(1,n);
energy_diff = zeros(1,n);
KOE = zeros(6,n);

for t = 1:n
    r_mag_vec(t) = sqrt( Y(1,t)^2 + Y(2,t)^2 + Y(3,t)^2 );
    v_mag_vec(t) = sqrt( Y(4,t)^2 + Y(5,t)^2 + Y(6,t)^2 );
    a_mag_vec(t) = -mu / r_mag_vec(t)^2;
    
    % calc change in specific energy over time
    energy_diff(t) = (v_mag_vec(t)^2 / 2 - mu /...
        r_mag_vec(t)) - energy0;
    
    % 6-by-n matrix of Keplarian components where rows correspond to
    % [a;e;i;RAAN;argp;tru] respectively
    KOE(:,t) = cart2kep_Whitbeck(mu, Y(:,t) );
end


% Plotting
figure
hold on

% Plotting magnitude of position over time
subplot(3,1,1)
plot(tspanhr,r_mag_vec)
title('Magnitude of Position of a Body in Orbit over Three Orbit Peiods')
xlabel('Time (hr)')
ylabel('Position (km)')

% Plotting magnitude of velocity over time
subplot(3,1,2)
plot(tspanhr,v_mag_vec)
title('Magnitude of Velocity of a Body in Orbit over Three Orbit Peiods')
xlabel('Time (hr)')
ylabel('Velocity (km/s)')

% Plotting magnitude of acceleration over time
subplot(3,1,3)
plot(tspanhr,a_mag_vec)
title('Magnitude of Acceleration of a Body in Orbit over Three Orbit Peiods')
xlabel('Time (hr)')
ylabel('Acceleration (km/s^2)')
hold off

%% #1 Part (b)
% Plotting Change in Specific Engergy of a Body in Orbit over Three Orbit Peiods
figure
hold on
plot(tspanhr,energy_diff)
title('Change in Specific Engergy of a Body in Orbit over Three Orbit Peiods')
xlabel('Time (hr)')
ylabel('Change in Specific Engergy (km^2/s)')
hold off

%% #1 Part (c)

%% Calc Keplarian orbital elements
% Semi-Major Axis (vector)
sma = KOE(1,:);
% Eccentricity (vector)
ecc = KOE(2,:);
% Inclination
inc = KOE(3,:);
% Right Ascension of Ascending Nodes
RAAN = KOE(4,:);
% Argument of Periapsis
argp = KOE(5,:);

%% Plotting Keplarian orbital elements
figure
hold on
set(0, 'DefaultAxesFontSize', 9)

% Semi-Major Axis
subplot(5,1,1)
plot(tspanhr,sma)
title('Semi-Major Axis (SMA) of a Body in Orbit over Three Orbit Peiods')
xlabel('Time (hr)')
ylabel('SMA (km)')

% Eccentricity
subplot(5,1,2)
plot(tspanhr,ecc)
title('Eccentricity of a Body in Orbit over Three Orbit Peiods')
xlabel('Time (hr)')
ylabel('Eccentricity')

% Inclination
subplot(5,1,3)
plot(tspanhr,inc)
title('Inclination of a Body in Orbit over Three Orbit Peiods')
xlabel('Time (hr)')
ylabel('Inclination (rad)')

% Right Ascension of Ascending Nodes
subplot(5,1,4)
plot(tspanhr,RAAN)
title('Right Ascension of Ascending Nodes (RAAN) of a Body in Orbit over Three Orbit Peiods')
xlabel('Time (hr)')
ylabel('RAAN (rad)')

% Argument of Periapsis
subplot(5,1,5)
plot(tspanhr,argp)
title('Argument of Periapsis (AoP) of a Body in Orbit over Three Orbit Peiods')
xlabel('Time (hr)')
ylabel('AoP (rad)')
hold off

%% #2: Factoring in J2 Perturbations
tspan2 = 0:20:86400; % seconds
tspan2hr = tspan2 / 3600; % hours

%% ODE45 Setup
% Anonymous function outputting 6x1 vector w/ 3D velocity & acceleration
% w/ input of time & 6x1 3D position/velocity vector.
vel_acc_J2 = @(t,r) [r(4:6);...
    -mu * r(1) / norm(r(1:3))^3 * (1 - J2 * 1.5 * (R/norm(r(1:3)))^2 * (5 * (r(3)/norm(r(1:3)))^2 - 1));...
    -mu * r(2) / norm(r(1:3))^3 * (1 - J2 * 1.5 * (R/norm(r(1:3)))^2 * (5 * (r(3)/norm(r(1:3)))^2 - 1));...
    -mu * r(3) / norm(r(1:3))^3 * (1 - J2 * 1.5 * (R/norm(r(1:3)))^2 * (5 * (r(3)/norm(r(1:3)))^2 - 3))];

% Calling ode45 function to output time vector T & position/velocity matrix
% at each time interval Y.
[T,Y2] = ode45(vel_acc_J2, tspan2, X0, odeoptions);

% Get Y matrix in form where columns contain pos/vel & each new column is
% another time t
Y2 = Y2';

% Initialize magnitudes of position, velocity, acceleration, energy change,
% & Kep elements vecs over time
n2 = length(tspan2);
r_mag_vec2 = zeros(1,n2);
v_mag_vec2 = zeros(1,n2);
a_mag_vec2 = zeros(1,n2);
energy_diff2 = zeros(1,n2);
energy2 = zeros(1,n2);
KOE2 = zeros(6,n2);

% Calc magnitudes of position, velocity & acceleration vecs over time
for t = 1:n2
    r_mag_vec2(t) = sqrt( Y2(1,t)^2 + Y2(2,t)^2 + Y2(3,t)^2 );
    v_mag_vec2(t) = sqrt( Y2(4,t)^2 + Y2(5,t)^2 + Y2(6,t)^2 );
    a_mag_vec2(t) = mu / r_mag_vec2(t)^2;
    
    % calc specific energy over time
    energy2(t) = (v_mag_vec2(t)^2 / 2) + (mu / r_mag_vec2(t)) ...
        - (mu / r_mag_vec2(t)) * (J2/2) * (R/r_mag_vec2(t))^2 ...
        * (3 * (Y2(3,t) / r_mag_vec2(t))^2 - 1);
        
    % calc change in specific energy over time factoring in J2 perturbations
    energy_diff2(t) = energy2(t) - energy2(1);
    
    % 6-by-n matrix of Keplarian components where rows correspond to
    % [a;e;i;RAAN;argp;tru] respectively
    KOE2(:,t) = cart2kep_Whitbeck(mu, Y2(:,t) );
end

% % Plotting
% figure
% hold on
% 
% % Plotting magnitude of position over time
% subplot(3,1,1)
% plot(tspan2hr,r_mag_vec2)
% title('Magnitude of Position of a Body in Orbit for One Earth Day Factoring in J2 Perturbations')
% xlabel('Time (hr)')
% ylabel('Position (km)')
% 
% % Plotting magnitude of velocity over time
% subplot(3,1,2)
% plot(tspan2hr,v_mag_vec2)
% title('Magnitude of Velocity of a Body in Orbit for One Earth Day Factoring in J2 Perturbations')
% xlabel('Time (hr)')
% ylabel('Velocity (km/s)')
% 
% % Plotting magnitude of acceleration over time
% subplot(3,1,3)
% plot(tspan2hr,a_mag_vec2)
% title('Magnitude of Acceleration of a Body in Orbit for One Earth Day Factoring in J2 Perturbations')
% xlabel('Time (hr)')
% ylabel('Acceleration (km/s^2)')
% hold off

%% #2 Part (b): Change in RAAN over time

%% Calc Keplarian orbital elements
% Semi-Major Axis (vector)
sma2 = KOE2(1,:);
% Eccentricity (vector)
ecc2 = KOE2(2,:);
% Inclination
inc2 = KOE2(3,:);
% Right Ascension of Ascending Nodes
RAAN2 = KOE2(4,:);
% Argument of Periapsis
argp2 = KOE2(5,:);

%% Plotting Keplarian orbital elements
figure
hold on

% Semi-Major Axis
subplot(5,1,1)
plot(tspan2hr,sma2)
title('Semi-Major Axis (SMA) of a Body in Orbit over Three Orbit Peiods')
xlabel('Time (hr)')
ylabel('SMA (km)')

% Eccentricity
subplot(5,1,2)
plot(tspan2hr,ecc2)
title('Eccentricity of a Body in Orbit over Three Orbit Peiods')
xlabel('Time (hr)')
ylabel('Eccentricity')

% Inclination
subplot(5,1,3)
plot(tspan2hr,inc2)
title('Inclination of a Body in Orbit over Three Orbit Peiods')
xlabel('Time (hr)')
ylabel('Inclination (rad)')

% Right Ascension of Ascending Nodes
subplot(5,1,4)
plot(tspan2hr,RAAN2)
title('Right Ascension of Ascending Nodes (RAAN) of a Body in Orbit over Three Orbit Peiods')
xlabel('Time (hr)')
ylabel('RAAN (rad)')

% Argument of Periapsis
subplot(5,1,5)
plot(tspan2hr,argp2)
title('Argument of Periapsis (AoP) of a Body in Orbit over Three Orbit Peiods')
xlabel('Time (hr)')
ylabel('AoP (rad)')
hold off

%% (2b). Calc change in RAAN over time dRAANdt
dRAANdt_eqn = - (3/2 * sqrt(mu) * J2 * R^2 / (1-ecc2(1)^2)^2 / sma2(1)^(7/2)) * cos(inc2(1)); % rad/s
% Convert to rad/day
dRAANdt_eqn = dRAANdt_eqn * 3600 * 24 % rad/day
% Find d(RAAN)/dt over 1 day from graph in 2a
RAAN_max = max(RAAN2);
RAAN_min = min(RAAN2);
dRAANdt_graph = RAAN_max - RAAN_min % rad/day
% Calc percent difference
dRAANdt_perc_diff = abs(dRAANdt_graph - dRAANdt_eqn) / dRAANdt_eqn * 100 % percent

%% Change in Specific Energy
% Calculated above
% Plotting
figure
hold on
plot(tspan2hr,energy_diff2)
title('Change in Specific Engergy of a Body in Orbit for One Earth Day Factoring in J2 Perturbations')
xlabel('Time (hr)')
ylabel('Change in Specific Engergy (km^2/s)')
hold off




