%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ASE 366L: Applied Orbital Mechanics
% HOMEWORK 2 MAIN FILE

close all; format long g; clc; clear all;
set(0, 'DefaultAxesFontSize', 10, 'defaultlinelinewidth',1)

%% Prob 1: Kep2Cart
% Convert given Keplarian Orbital Elements (KOE) to Cartesian elements
mu_earth = 398600.4415; % standard gravitation parameter of Earth (km^3/s^2)
% xKOE follows format: <SMA,Ecc,Inc,RAAN,ArgPeri,t_p>
KOE0 = [26000, .72, 75*pi/180, 90*pi/180, -90*pi/180, 0];

disp('Problem 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
[Cart0] = kep2cart_Whitbeck( mu_earth, KOE0 )

%% Prob 2:
sma = 26000; % km
ecc = 0.72;
inc = 75*pi/180; % rads
raan = 90*pi/180; % rads
arg = -90*pi/180; % rads

t_p = 0; % Time since last Periapse passage (second)
t_f = 3*86400; % seconds
t_vec = (0:20:t_f); % seconds

% Calls function that propagates the KOEs via true anomaly
[KOE_mat] = prop_KOE(sma,ecc,inc,raan,arg,t_p,t_vec,mu_earth);

% Displays final true anomaly in degrees
disp('Problem 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
tru_final = KOE_mat(length(t_vec),6)*180/pi

% Displays final mean anomaly in degrees
Me_final = KOE_mat(length(t_vec),7)*180/pi

%% Prob 3:
% Convert propagated KOEs to Cartesian elements at each time

% Get necessary # of iterations of orbit propagation to use in Kep2Cart
[m,n] = size(KOE_mat);

% Initializing matrix of Cartesian elements
Cart_mat = zeros(m,6);

% Calculating matrix of Cartesian elements
for ii = 1:m
   Cart_mat(ii,:) = kep2cart_Whitbeck(mu_earth,KOE_mat(ii,1:6));
end

% Displays final position vector (km)
disp('Problem 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
position_final = Cart_mat(m,1:3)

% Displays final velocity vector (km/s)
velocity_final = Cart_mat(m,4:6)


%% Prob 4: Propagate Cart0 numerically

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Givens %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
vel_acc = @(t,r) [r(4:6); -mu_earth .* r(1:3) ./ (norm(r(1:3))^3)];

% Time span every 60 seconds starting w/ initial time t0 & ending w/ final
% time tf.
t0 = 0; % seconds
tf = 3*86400; % seconds
tspan = (t0:20:tf); % seconds

% ODE options as given in lecture
odeoptions = odeset('RelTol',1e-12,'AbsTol',1e-20);

% Calling ode45 function to output time vector T & position/velocity matrix
% at each time interval Y.
[~,Y] = ode45(vel_acc, tspan', Cart0, odeoptions);

% Get Y matrix in form where columns contain pos/vel & each new column is
% another time t
Y = Y';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculating Specific Energies %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize magnitudes of position, velocity & acceleration vecs over time
n = length(tspan);
r_mag_vec = zeros(1,n);
v_mag_vec = zeros(1,n);
a_mag_vec = zeros(1,n);
energy_diff = zeros(1,n);
KOE = zeros(6,n);

% Calc magnitudes of position, velocity & acceleration vecs over time
for t = 1:n
    r_mag_vec(t) = sqrt( Y(1,t)^2 + Y(2,t)^2 + Y(3,t)^2 );
    v_mag_vec(t) = sqrt( Y(4,t)^2 + Y(5,t)^2 + Y(6,t)^2 );
    a_mag_vec(t) = -mu_earth / r_mag_vec(t)^2;
    
    % calc change in specific energy over time
    energy_diff(t) = (v_mag_vec(t)^2 / 2 - mu_earth /...
        r_mag_vec(t)) - energy0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting Change in Specific Engergy of a Body in Orbit %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
hold on
grid on
plot(tspan/3600,energy_diff)
title('Change in Specific Engergy of a Body in Orbit')
xlabel('Time (hours)')
ylabel('Change in Specific Engergy (km^2/s^2)')
hold off

%% Prob 5: Error between Final Propagations of ODE45 & True Anomaly

% ODE45 final position vector
ODE_pos_final = Y(1:3,length(tspan))'; % km
% ODE45 final velocity vector
ODE_vel_final = Y(4:6,length(tspan))'; % km/s

disp('Problem 5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
% Error between final positions
Error_pos_final = ODE_pos_final - position_final % km
% Error between final velocities
Error_vel_final = ODE_vel_final - velocity_final % km/s

%% Prob 6: Plotting Error btwn entire Propagations for ODE45 & True Anomaly

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Position Error
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 3D position for True Anomaly propagation
Tru_pos_mat = Cart_mat(:,1:3);
% 3D position for ODE45 propagation
ODE_pos_mat = Y(1:3,:)';
% Error for 3D position
error_pos_mat = ODE_pos_mat - Tru_pos_mat;

% Separating the position error dimensions
error_x_pos = error_pos_mat(:,1);
error_y_pos = error_pos_mat(:,2);
error_z_pos = error_pos_mat(:,3);

% Plotting the position error dimensions
figure
hold on

subplot(3,1,1)
plot(tspan/3600,error_x_pos)
title('X-dimension Position Error between entire Propagations for ODE45 & True Anomaly')
xlabel('Time (hours)')
ylabel('X-Position Error (km)')
grid on

subplot(3,1,2)
plot(tspan/3600,error_y_pos)
title('Y-dimension Position Error between entire Propagations for ODE45 & True Anomaly')
xlabel('Time (hours)')
ylabel('Y-Position Error (km)')
grid on

subplot(3,1,3)
plot(tspan/3600,error_z_pos)
title('Z-dimension Position Error between entire Propagations for ODE45 & True Anomaly')
xlabel('Time (hours)')
ylabel('Z-Position Error (km)')
grid on

hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Velocity Error
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 3D Velocity for True Anomaly propagation
Tru_vel_mat = Cart_mat(:,4:6);
% 3D Velocity for ODE45 propagation
ODE_vel_mat = Y(4:6,:)';
% Error for 3D Velocity
error_vel_mat = ODE_vel_mat - Tru_vel_mat;

% Separating the Velocity error dimensions
error_x_vel = error_vel_mat(:,1);
error_y_vel = error_vel_mat(:,2);
error_z_vel = error_vel_mat(:,3);

% Plotting the Velocity error dimensions
figure
hold on

subplot(3,1,1)
plot(tspan/3600,error_x_vel)
title('X-dimension Velocity Error between entire Propagations for ODE45 & True Anomaly')
xlabel('Time (hours)')
ylabel('X-Velocity Error (km/s)')
grid on

subplot(3,1,2)
plot(tspan/3600,error_y_vel)
title('Y-dimension Velocity Error between entire Propagations for ODE45 & True Anomaly')
xlabel('Time (hours)')
ylabel('Y-Velocity Error (km/s)')
grid on

subplot(3,1,3)
plot(tspan/3600,error_z_vel)
title('Z-dimension Velocity Error between entire Propagations for ODE45 & True Anomaly')
xlabel('Time (hours)')
ylabel('Z-Velocity Error (km/s)')
grid on

hold off