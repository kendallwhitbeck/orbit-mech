%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ASE 366L: Applied Orbital Mechanics
% HOMEWORK 7 MAIN FILE
tic
close all; format long g; clc; clear all;
% set(0, 'DefaultAxesFontSize', 10, 'defaultlinelinewidth',1)

% Problem 3
mu_earth = 398600.4415; % km^3/s^2
xKOE = [7000, .01, pi/4, 0, pi/2, -pi/2]; 

Cart0 = kep2cart_Whitbeck(mu_earth,xKOE);
r_sc_vec = Cart0(1:3); % km
v_sc_vec = Cart0(4:6); % km/s
r_sc = norm(r_sc_vec); % magnitude of position vector (km)
v_sc = norm(v_sc_vec); % magnitude of position vector (km/s)

% Part (b.)
R_earth = 6378.1363; % radius of earth (km)
J2 = 0.0010826267; % J2 perturbation

r_i = r_sc_vec(1); % i-component of position vector (km)
r_j = r_sc_vec(2); % j-component of position vector (km)
r_k = r_sc_vec(3); % k-component of position vector (km)

% J2
a_J2_i = 3*mu_earth*J2*R_earth^2 / (2*r_sc^5) * r_i * (5 * (r_k/r_sc)^2 - 1);
a_J2_j = 3*mu_earth*J2*R_earth^2 / (2*r_sc^5) * r_j * (5 * (r_k/r_sc)^2 - 1);
a_J2_k = 3*mu_earth*J2*R_earth^2 / (2*r_sc^5) * r_k * (5 * (r_k/r_sc)^2 - 3);

a_J2_vec = [a_J2_i, a_J2_j, a_J2_k];

disp('The answer to #3, Part (b.):')
a_J2 = norm(a_J2_vec) % km/s^2

% Part (c)
C_D = 2.0;
A_3c = 1.5; % m^2
m_3c = 1000; % kg
A2m_3c = A_3c / m_3c; % m^2/kg

a_D_vec_3c = accel_drag(Cart0,C_D,A2m_3c);
disp('The answer to #3, Part (c.):')
a_D_3c = norm(a_D_vec_3c)

%% Problem 4
h_range = (0:1:1000); % Range of altitudes (km)
r_range = h_range + R_earth; % Range of orbit radii (km)
rho_range = zeros(1,length(r_range));

for ii = 1:length(r_range)
    rho_range(ii) = atmo_density(r_range(ii));
end

figure
semilogx(rho_range,h_range)
grid on
title('U.S. Standard Atmospheric Density for Varying Altitude')
xlabel('Atmospheric Density (kg/km^3)')
ylabel('Altitude (km)')

%% Problem 5
% Givens
SMA = 7300.0; % km
ecc = 0.01;
inc = 45*pi/180; % rad
RAAN = 0; % rad
argp = 0; % rad
tru = 0; % rad
KOE0 = [SMA,ecc,inc,RAAN,argp,tru];

C_D = 2.0; % Coefficient of Drag
A2m = 0.01; % Area-to-Mass Ratio (m^2/kg)

sec2day = 1/3600/24; % Seconds-to-Days Conversion Factor (days/sec)
tspan = (0 : 0.5/sec2day : 100/sec2day); % 10 days every 0.5 days (in sec)
[m,n] = size(tspan);

% Convert initial KOEs into Cartesian position/velocity
[Cart0] = kep2cart_Whitbeck( mu_earth, KOE0 );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ODE45 Setup %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Anonymous function outputting 6x1 vector w/ 3D velocity & acceleration
% w/ input of time & 6x1 3D position/velocity vector.
vel_acc = @(t,Cart) [Cart(4:6);...
    -mu_earth .* Cart(1:3) ./ (norm(Cart(1:3))^3)...
    + accel_drag(Cart,C_D,A2m)];

% ODE options as given in lecture
odeoptions = odeset('RelTol',1e-12,'AbsTol',1e-30);

% Calling ode45 function to output time vector T & position/velocity matrix
% at each time interval Y.
[~,Y] = ode45(vel_acc, tspan, Cart0, odeoptions);

% Get Y matrix in form where columns contain pos/vel & each new column is
% another time t
Y = Y';
r_t_mat = Y(1:3,:);
v_t_mat = Y(4:6,:);

% Initial Position vector
r0 = Cart0(1:3); % km
% Initial Velocity vector
v0 = Cart0(4:6); % km/s
% Initial Energy scalar using Vis-Viva
energy0 = dot(v0,v0) / 2 - mu_earth / norm(r0);

% Calculation of values to be plotted
radius = zeros(1,length(tspan));
speed = zeros(1,length(tspan));
energy_diff = zeros(1,length(tspan));
accel_vec = zeros(3,length(tspan));
accel_mag = zeros(1,length(tspan));

for ii = 1:length(tspan)
    radius(ii) = norm(r_t_mat(:,ii));
    speed(ii) = norm(v_t_mat(:,ii));
        
    % calc change in specific energy at each time iteration
    energy_diff(ii) = (speed(ii)^2 / 2 - mu_earth / radius(ii)) - energy0;
    
    % Calc Keplarian Elements at each Time iteration
    KOE_range(:,ii) = cart2kep_Whitbeck(mu_earth, Y(:,ii));
    SMA_range(ii) = KOE_range(1,ii);
    ecc_range(ii) = KOE_range(2,ii);

    % Calculate Radius of Periapse and Apoapse
    r_p(ii) = SMA_range(ii) * (1 - ecc_range(ii));
    r_a(ii) = SMA_range(ii) * (1 + ecc_range(ii));

end

%% Problem 5, Part (a)
figure
hold on
grid on
plot(tspan.*sec2day,energy_diff)
title('Change in Specific Engergy over Time')
xlabel('Time (days)')
ylabel('Change in Specific Engergy (km^2/s^2)')
hold off

%% Problem 5, Part (b)
figure
hold on

subplot(3,1,1)
plot(tspan.*sec2day,radius)
title('Orbit Radius over Time')
xlabel('Time (days)')
ylabel('Orbit Radius (km)')
grid on

subplot(3,1,2)
plot(tspan.*sec2day,r_p)
title('Radius of Periapse over Time')
xlabel('Time (days)')
ylabel('Radius of Periapse (km)')
grid on

subplot(3,1,3)
plot(tspan.*sec2day,r_a)
title('Radius of Apoapse over Time')
xlabel('Time (days)')
ylabel('Radius of Apoapse (km)')
grid on

hold off

%% Problem 5, Part (c)
figure
hold on

subplot(2,1,1)
plot(tspan.*sec2day,SMA_range)
title('Semi-Major Axis over Time')
xlabel('Time (days)')
ylabel('Semi-Major Axis (km)')
grid on

subplot(2,1,2)
plot(tspan.*sec2day,ecc_range)
title('Eccentricity over Time')
xlabel('Time (days)')
ylabel('Eccentricity')
grid on

hold off

time_to_run = toc




