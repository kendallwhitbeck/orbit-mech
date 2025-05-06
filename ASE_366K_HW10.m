%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ASE 366K: Spacecraft Dynamics
% HOMEWORK 10 MAIN FILE
tic
close all; format shortg; clc; % clear all;
set(0, 'DefaultAxesFontSize', 16, 'defaultlinelinewidth',1)

%% Prob 1: Generate Groundtracks
mu = 398600; % standard gravitation parameter of Earth (km^3/s^2)
theta_G0 = 0; % initial Greenwich sidereal time (rad)
f = 3.35e-3; % oblateness
angvel_planet = 2*pi / 86164; % angular velocity of Earth (rad/s)

%% Part (a)
% Givens
SMA1 = 10000; % Semi-Major Axis (km)

%% Convert given Keplarian Orbital Elements (KOE) to Cartesian elements

xKOEa = [10000, .5, 85*pi/180, 0, -55*pi/180, 0];
[xC0] = kep2cart_Whitbeck( mu, xKOEa );

%% Propagate xC0

% Anonymous function outputting 6x1 vector w/ 3D velocity & acceleration
% w/ input of time & 6x1 3D position/velocity vector.
vel_acc = @(t,r) [r(4:6); -mu .* r(1:3) ./ (norm(r(1:3))^3)];

% Time span every 60 seconds starting w/ initial time t0 & ending w/ final
% time tf.
t0 = 0; % seconds
tf = 29856; % seconds (=.3465 days)
tspan = (t0:60:tf); % seconds

% ODE options as given in lecture
odeoptions = odeset('RelTol',1e-15,'AbsTol',1e-20);

% Calling ode45 function to output time vector T & position/velocity matrix
% at each time interval Y.
[Tvec,Y] = ode45(vel_acc, tspan, xC0, odeoptions);

% Pulling position data values (x,y,z) from Y matrix
rmat_ECI = Y(:,1:3);

%% Convert Cartesian Propagation to Geodetic Latitude & Longitude
% These are in RADIANS!!!
[lat_gd_rad, long_rad] = cart2lat(t0, theta_G0, angvel_planet, rmat_ECI, Tvec, f);

% Convert to degrees
lat_gd_deg = lat_gd_rad*180/pi;
long_deg = long_rad*180/pi;

%% Plotting the Groundtracks

% Loading the earth coastline plot data
load earth_coastline.mat

% Plotting the Coastline and Ground Tracks
figure
hold on
grid on
plot(earth_coastline(:,1),earth_coastline(:,2),'k')
plot(-180:180,zeros(361),'--','Color','k','LineWidth',.25)
plot(long_deg,lat_gd_deg,'o','MarkerEdgeColor','r','LineWidth',1.5)
axis equal
xlim([-180,180])
ylim([-90,90])
xlabel('Longitude (deg)')
ylabel('Geodetic Latitude (deg)')
set(gca,'Xtick',-150:50:150)
set(gca,'Ytick',-80:20:80)
hold off

% Calculating three periods using phase shift from groundtracks
phase_deg_a = abs(-48.45 - -7.095);
period_groundtracks_a = phase_deg_a * pi / 180 / ( 2*pi/86164 )

% Formulaically calculating three periods
period_formula_a = 2*pi*sqrt(SMA1^3/mu) % seconds

% Calc difference between the two
period_diff_percent_a = abs(period_groundtracks_a - period_formula_a)/period_formula_a*100

%% Part (b)
% Givens
SMA2 = 50000; % Semi-Major Axis (km)

%% Convert given Keplarian Orbital Elements (KOE) to Cartesian elements

xKOEa = [50000, 0, 45*pi/180, 100*pi/180, 90*pi/180, 0];
[xC0] = kep2cart_Whitbeck( mu, xKOEa );

%% Propagate xC0

% Anonymous function outputting 6x1 vector w/ 3D velocity & acceleration
% w/ input of time & 6x1 3D position/velocity vector.
vel_acc = @(t,r) [r(4:6); -mu .* r(1:3) ./ (norm(r(1:3))^3)];

% Time span every 60 seconds starting w/ initial time t0 & ending w/ final
% time tf.
t0 = 0; % seconds
tf = 222530; % seconds
tspan = (t0:60:tf); % seconds

% ODE options as given in lecture
odeoptions = odeset('RelTol',1e-13,'AbsTol',1e-20);

% Calling ode45 function to output time vector T & position/velocity matrix
% at each time interval Y.
[Tvec,Y] = ode45(vel_acc, tspan, xC0, odeoptions);

% Pulling position data values (x,y,z) from Y matrix
rmat_ECI = Y(:,1:3);

%% Convert Cartesian Propagation to Geodetic Latitude & Longitude
% These are in RADIANS!!!
[lat_gd_rad, long_rad] = cart2lat(t0, theta_G0, angvel_planet, rmat_ECI, Tvec, f);

% Convert to degrees
lat_gd_deg = lat_gd_rad*180/pi;
long_deg = long_rad*180/pi;

%% Plotting the Groundtracks

% Loading the earth coastline plot data
load earth_coastline.mat

% plot everything
figure
hold on
grid on
plot(earth_coastline(:,1),earth_coastline(:,2),'k')
plot(-180:180,zeros(361),'--','Color','k','LineWidth',.25)
plot(long_deg,lat_gd_deg,'o','MarkerEdgeColor','r')
axis equal
xlim([-180,180])
ylim([-90,90])
xlabel('Longitude (deg)')
ylabel('Geodetic Latitude (deg)')
set(gca,'Xtick',-150:50:150)
set(gca,'Ytick',-80:20:80)
hold off

% Calculating two periods using phase shift from groundtracks
phase_deg_b = abs(-4.831 - 100) + 360
period_groundtracks_b = phase_deg_b * pi / 180 / ( 2*pi/86164 )

% Formulaically calculating three periods
period_formula_b = 2*pi*sqrt(SMA2^3/mu) % seconds

% Calc difference between the two periods
period_diff_percent_b = abs(period_groundtracks_b - period_formula_b)/period_formula_b*100

%% Part (c)
% Givens
M = 5; N = 6;

% Calc SMA from repeating orbit ground track
SMA3 = (mu * ( M / N / angvel_planet )^2 )^(1/3) % Semi-Major Axis (km)

%% Convert given Keplarian Orbital Elements (KOE) to Cartesian elements

xKOEa = [SMA3, 0, 55*pi/180, 200*pi/180, 0*pi/180, 0];
[xC0] = kep2cart_Whitbeck( mu, xKOEa );

%% Propagate xC0

% Anonymous function outputting 6x1 vector w/ 3D velocity & acceleration
% w/ input of time & 6x1 3D position/velocity vector.
vel_acc = @(t,r) [r(4:6); -mu .* r(1:3) ./ (norm(r(1:3))^3)];

% Time span every 60 seconds starting w/ initial time t0 & ending w/ final
% time tf.
t0 = 0; % seconds
tf = 5.0262e+05; % seconds
tspan = (t0:60:tf); % seconds

% ODE options as given in lecture
odeoptions = odeset('RelTol',1e-13,'AbsTol',1e-20);

% Calling ode45 function to output time vector T & position/velocity matrix
% at each time interval Y.
[Tvec,Y] = ode45(vel_acc, tspan, xC0, odeoptions);

% Pulling position data values (x,y,z) from Y matrix
rmat_ECI = Y(:,1:3);

%% Convert Cartesian Propagation to Geodetic Latitude & Longitude
% These are in RADIANS!!!
[lat_gd_rad, long_rad] = cart2lat(t0, theta_G0, angvel_planet, rmat_ECI, Tvec, f);

% Convert to degrees
lat_gd_deg = lat_gd_rad*180/pi;
long_deg = long_rad*180/pi;

%% Plotting the Groundtracks

% Loading the earth coastline plot data
load earth_coastline.mat

% plot everything
figure
hold on
grid on
plot(earth_coastline(:,1),earth_coastline(:,2),'k')
plot(-180:180,zeros(361),'--','Color','k','LineWidth',.25)
plot(long_deg,lat_gd_deg,'o','MarkerEdgeColor','r')
axis equal
xlim([-180,180])
ylim([-90,90])
xlabel('Longitude (deg)')
ylabel('Geodetic Latitude (deg)')
set(gca,'Xtick',-150:50:150)
set(gca,'Ytick',-80:20:80)
hold off

% Formulaically calculating period to be used in timespan
period_formula_c = 2*pi*sqrt(SMA3^3/mu); % seconds

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
toc