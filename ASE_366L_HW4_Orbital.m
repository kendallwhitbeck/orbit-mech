%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ASE 366L: Applied Orbital Mechanics
% HOMEWORK 4 MAIN FILE

close all; format long g; clc; clear all;
% set(0, 'DefaultAxesFontSize', 10, 'defaultlinelinewidth',1)

% Given:
dUT1 = 0; % UTC to UT1 conversion factor

%% Prob 4
t_test = 2451545.0; % days (J2000.0, UTC)

% Part (a.)
JD_UTC_0 = 2457793.5; % Initial time in days (J2000.0, UTC)

% Conversion factor for Astronimcal Units to kilometers
AU2km = 149597870.7;

% USing sun function to get sun's position in J2000.0 coordinates
r_sun_J2000_vec = sun(JD_UTC_0,dUT1) % AU

% Part (b.)
t_0 = JD_UTC_0; % days
t_f = t_0 + 365.25; % days
tspan = linspace(t_0,t_f,365); % days

figure
hold on
for ii = 1:length(tspan)
    
    r_sun_J2000_vec(:,ii) = sun(tspan(ii),dUT1); % AU
    
end

subplot(3,1,1)
plot(1:length(tspan),r_sun_J2000_vec(1,:))
title('X-Component of Sun''s Position in J2000.0 Frame over One Year')
xlabel('Time (days)')
ylabel('Sun''s X-Position (AU)')

subplot(3,1,2)
plot(1:length(tspan),r_sun_J2000_vec(2,:))
title('Y-Component of Sun''s Position in J2000.0 Frame over One Year')
xlabel('Time (days)')
ylabel('Sun''s Y-Position (AU)')

subplot(3,1,3)
plot(1:length(tspan),r_sun_J2000_vec(3,:))
title('Z-Component of Sun''s Position in J2000.0 Frame over One Year')
xlabel('Time (days)')
ylabel('Sun''s Z-Position (AU)')

hold off

%% Prob 5

% Givens
mu_earth = 3998600.4415; % km^3/s^2
mu_sun = 1.3271244e11; % km^3/s^2

t_0 = 245154.; % start time in days (J2000.0 UTC)
t_f = t_0 + 365.25; % days
sec2day = 1/3600/24; % days per second

tspan_days = linspace(0, 365.25, 365*2); % days
tspan_secs = tspan_days .* 24 .* 3600; % seconds

% Initial KOEs
sma = 500000.000; % km
ecc = 0.1;
inc = 50*pi/180; % rads
raan = 350*pi/180; % rads
arg = 90*pi/180; % rads
tru = 0; % rads
KOE0 = [sma,ecc,inc,raan,arg,tru];

% Convert initial KOEs into Cartesian position/velocity
[Cart0] = kep2cart_Whitbeck( mu_earth, KOE0 )

% Inertial, geocentric, ijk frame at time t0 = 0
% Initial Position vector
r0 = Cart0(1:3); % km
% Initial Velocity vector
v0 = Cart0(4:6); % km/s
% Initial Energy scalar using Vis-Viva
energy0 = dot(v0,v0) / 2 - mu_earth / norm(r0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ODE45 Setup %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Anonymous function outputting 6x1 vector w/ 3D velocity & acceleration
% w/ input of time & 6x1 3D position/velocity vector.
vel_acc = @(t,r) [r(4:6);
    -mu_earth .* r(1:3) ./ (norm(r(1:3))^3) + mu_sun...
    .* ( (r(1:3)+AU2km*sun(t*sec2day+t_0,0))...
    ./(norm(r(1:3)+AU2km*sun(t*sec2day+t_0,0))^3)...
    - AU2km*sun(t*sec2day+t_0,0)./(norm(AU2km*sun(t*sec2day+t_0,0))^3) )];

% ODE options as given in lecture
odeoptions = odeset('RelTol',1e-18,'AbsTol',1e-30);

% Calling ode45 function to output time vector T & position/velocity matrix
% at each time interval Y.
[~,Y] = ode45(vel_acc, tspan_secs', Cart0, odeoptions);

% Get Y matrix in form where columns contain pos/vel & each new column is
% another time t
Y = Y';

% Convert each Cartesian element of Y matrix to KOE
for ii = 1:length(tspan_secs)
    
    KOE_prop(:,ii) = cart2kep_Whitbeck(mu_earth,Y(:,ii));

end
%%
figure
hold on

subplot(3,2,1)
plot(tspan_days,KOE_prop(1,:))
title('Semi-Major Axis (SMA) vs. Time')
xlabel('Time (days)')
ylabel('SMA (km)')
xlim([0 367])

subplot(3,2,2)
plot(tspan_days,KOE_prop(2,:))
title('Eccentricity vs. Time')
xlabel('Time (days)')
ylabel('Eccentricity')
xlim([0 367])

subplot(3,2,3)
plot(tspan_days,KOE_prop(3,:))
title('Inclination vs. Time')
xlabel('Time (days)')
ylabel('Inclination (rads)')
xlim([0 367])

subplot(3,2,4)
plot(tspan_days,KOE_prop(4,:))
title('Right Ascension of the Ascending Node (RAAN) vs. Time')
xlabel('Time (days)')
ylabel('RAAN (rads)')
xlim([0 367])

subplot(3,2,5)
plot(tspan_days,KOE_prop(5,:))
title('Argument of Periapsis (ArgP) vs. Time')
xlabel('Time (days)')
ylabel('ArgP (rads)')
xlim([0 367])

subplot(3,2,6)
plot(tspan_days,KOE_prop(6,:))
title('True Anomaly (Tru) vs. Time')
xlabel('Time (days)')
ylabel('Tru (rad)')
xlim([0 367])

hold off