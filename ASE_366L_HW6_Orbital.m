%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ASE 366L: Applied Orbital Mechanics
% HOMEWORK 6 MAIN FILE

close all; format long g; clc; clear all;
% set(0, 'DefaultAxesFontSize', 10, 'defaultlinelinewidth',1)

%% Problem 5
r_vec = [.009, 4595.737, 4595.731]; % km

disp('Problem 5:')
disp('The perturbing acceleration resulting from J2 & J3 in km/s^2 is:')
[a_p_vec] = accel_J2_J3(r_vec)

%% Problem 6, Part (a.)

% Givens
mu_earth = 398600.4415; % gravitational parameter (km^3/s^2)
R_earth = 6378.1363; % radius of earth (km)

SMA = 7000; % km
ecc = .01;
inc = pi/4; % rad
RAAN = 0; % rad
argp = -pi/2; % rad
tru = 2*pi/3; % rad
KOE0 = [SMA,ecc,inc,RAAN,argp,tru];

tspan = linspace(0,24*3600,100); % 24 hours represented in seconds
[~,n] = size(tspan);

% Convert initial KOEs into Cartesian position/velocity
[Cart0] = kep2cart_Whitbeck( mu_earth, KOE0 );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ODE45 Setup %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Anonymous function outputting 6x1 vector w/ 3D velocity & acceleration
% w/ input of time & 6x1 3D position/velocity vector.
vel_acc = @(t,Cart) [Cart(4:6);...
    -mu_earth .* Cart(1:3) ./ (norm(Cart(1:3))^3)...
    + accel_J2_J3(Cart(1:3))'];

% ODE options as given in lecture
odeoptions = odeset('RelTol',1e-18,'AbsTol',1e-30);

% Calling ode45 function to output time vector T & position/velocity matrix
% at each time interval Y.
[~,Y] = ode45(vel_acc, tspan', Cart0, odeoptions);

% Get Y matrix in form where columns contain pos/vel & each new column is
% another time t
Y = Y';

disp('Problem 6, Part (a):')
disp('The final position vector after 24 hours in kilometers is:')
r_final_vec = Y(1:3,n)

%% Problem 6, Part (b.)
% Givens
SMA = 6800.0; % km
ecc = 0.03;
inc = 1.309; % rad
RAAN = 3.49; % rad
argp = 5.24; % rad
tru = 0; % rad
KOE0 = [SMA,ecc,inc,RAAN,argp,tru];

tspan = (0:60:24*3600); % 24 hours represented in seconds every 60 seconds
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
    + accel_J2_J3(Cart(1:3))'];

% ODE options as given in lecture
odeoptions = odeset('RelTol',1e-18,'AbsTol',1e-30);

% Calling ode45 function to output time vector T & position/velocity matrix
% at each time interval Y.
[~,Y] = ode45(vel_acc, tspan', Cart0, odeoptions);

% Get Y matrix in form where columns contain pos/vel & each new column is
% another time t
Y = Y';
r_t_mat = Y(1:3,:);
v_t_mat = Y(4:6,:);

disp('Problem 6, Part (b)')
disp('The final position vector after 24 hours in kilometers is:')
r_final_vec = r_t_mat(:,n)

% Calculation of values to be plotted
for ii = 1:length(tspan)
    radius(ii) = norm(r_t_mat(:,ii));
    speed(ii) = norm(v_t_mat(:,ii));
    
    accel_vec(1:3,ii) = -mu_earth .* r_t_mat(:,ii) ./ ...
        radius(ii)^3 + accel_J2_J3(r_t_mat(:,ii))';
    
    accel_mag(ii) = norm(accel_vec(:,ii));
end

% Plotting
figure
hold on

subplot(3,1,1)
plot(tspan/3600,radius)
title('Radius of Satellite over 24 Hour Timespan')
xlabel('Time (hours)')
ylabel('Radius (km)')
xlim([0 24])
grid on

subplot(3,1,2)
plot(tspan/3600,speed)
title('Speed of Satellite over 24 Hour Timespan')
xlabel('Time (hours)')
ylabel('Speed (km/s)')
xlim([0 24])
grid on

subplot(3,1,3)
plot(tspan/3600,accel_mag)
title('Magnitude of Acceleration of Satellite over 24 Hour Timespan')
xlabel('Time (hours)')
ylabel('Magnitude of Acceleration (km/s^2)')
xlim([0 24])
grid on

