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
% Given KOEs for target
sma = 7000; % km
ecc = 0;
inc = 75*pi/180; % rad
RAAN = 0; % rad
argp = 0; % rad
tru = 0; % rad

% Vector of initial Keplarian Orbital Elements for target
xKOE = [sma, ecc, inc, RAAN, argp, tru];

% Convert initial KOEs into Cartesian position/velocity
Cart0_targ_ijk = kep2cart_Whitbeck(mu_earth,xKOE);

% Inertial (ijk frame) position & velocity of the target
r_targ_ijk = Cart0_targ_ijk(1:3); % km
v_targ_ijk = Cart0_targ_ijk(4:6); % km/s

% Initial Relative (RSW frame) Position & Velocity of Chaser
rho_0_chas_RSW = [25; 0; 0] ./ 1000; % km
zeta_0_chas_RSW = [-1; 0; 0] ./ 1000; % km/s

% Transform Relative position from RSW to ijk
% Calculating necessary axes
R_unit = [r_targ_ijk / norm(r_targ_ijk)];
W_unit = [cross(r_targ_ijk,v_targ_ijk) ./ (norm(r_targ_ijk) * norm(v_targ_ijk))];
S_unit = cross(W_unit,R_unit);

% Defining Transformation Matrix for RSW -> ijk
T_RSW_ijk = [R_unit, S_unit, W_unit];

% Transforming relative position & velocity from RSW to ijk
rho_0_chas_ijk = T_RSW_ijk * rho_0_chas_RSW;
zeta_0_chas_ijk = T_RSW_ijk * zeta_0_chas_RSW;

% Calc chaser position & velocity in inertial frame
r_chas_ijk = rho_0_chas_ijk + r_targ_ijk; % km
v_chas_ijk = zeta_0_chas_ijk + v_targ_ijk; % km/s

fprintf('-------------- PROBLEM 5, PART (b) -------------- \n\n')
fprintf('The Inertial Position of the Chaser in km is:\n')
disp(r_chas_ijk)
fprintf('The Inertial Velocity of the Chaser in km/s is:\n')
disp(v_chas_ijk)

%% Problem 5, Part (c.)

% Calculate 3 Orbit Periods, T3, & create matching timespan
T3 = 3 * ceil(2 * pi * sqrt(sma^3 / mu_earth)); % sec
tspan = (0:20:T3); % 3 orbit periods (sec)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ODE45 Setup %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Anonymous function outputting 6x1 vector w/ 3D velocity & acceleration
% w/ input of time & 6x1 3D position/velocity vector.
vel_acc = @(t,Cart) [Cart(4:6);...
    -mu_earth .* Cart(1:3) ./ (norm(Cart(1:3))^3)];

% ODE options as given in lecture
odeoptions = odeset('RelTol',2.22045e-14,'AbsTol',1e-30);

% Calling ode45 function to output time vector T & position/velocity matrix
% at each time interval Y.
[~,Y_targ_ijk] = ode45(vel_acc, tspan, Cart0_targ_ijk, odeoptions);

% Get Y matrix in form where columns contain pos/vel & each new column is
% another time t, then break it into pos matrix & vel matrix
Y_targ_ijk = Y_targ_ijk';

% Define inertial Cartesian state for chaser
Cart0_chas_ijk = [r_chas_ijk; v_chas_ijk]; % km & km/s

% Calling ode45 function to output time vector T & position/velocity matrix
% at each time interval Y.
[~,Y_chas_ijk] = ode45(vel_acc, tspan, Cart0_chas_ijk, odeoptions);

% Get Y matrix in form where columns contain pos/vel & each new column is
% another time t, then break it into pos matrix & vel matrix
Y_chas_ijk = Y_chas_ijk';

% Position vectors over time for target & chaser
r_targ_mat_ijk = Y_targ_ijk(1:3,:);
r_chas_mat_ijk = Y_chas_ijk(1:3,:);

% Velocity vectors over time for target & chaser
v_targ_mat_ijk = Y_targ_ijk(4:6,:);
v_chas_mat_ijk = Y_chas_ijk(4:6,:);

% Relative position of Chaser to target over time in ijk frame
rho_chas_mat_ijk = r_chas_mat_ijk - r_targ_mat_ijk;

rho_chas_mat_RSW = zeros(3,length(tspan));
% needs new T_RSW_ijk for each r_targ_ijk %%%%%%%
for ii = 1:length(tspan)
    
    % Calculating necessary axes
    R_unit = [r_targ_mat_ijk(:,ii) / norm(r_targ_mat_ijk(:,ii))];
    W_unit = [cross(r_targ_mat_ijk(:,ii),v_targ_mat_ijk(:,ii))...
        ./ (norm(r_targ_mat_ijk(:,ii)) * norm(v_targ_mat_ijk(:,ii)))];
    S_unit = cross(W_unit,R_unit);
    
    % Defining Transformation Matrix for RSW -> ijk
    T_RSW_ijk = [R_unit, S_unit, W_unit];
    % Transposing to get transformation matrix for ijk -> RSW
    T_ijk_RSW = T_RSW_ijk';
    
    % Transforming relative position of chaser from ijk to RSW frame
    rho_chas_mat_RSW(:,ii) = T_ijk_RSW * rho_chas_mat_ijk(:,ii);
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% R-component over time
rho_R = rho_chas_mat_RSW(1,:); % km
% S-component over time
rho_S = rho_chas_mat_RSW(2,:); % km

% Plot rho_R vs. rho_S
figure
hold on
plot(rho_S*1000,rho_R*1000)
grid on
title('5c) R-component vs. S-component of Relative Position for a Chaser & Target Over Three Orbit Periods')
xlabel('\rho_S (m)')
ylabel('\rho_R (m)')
hold off

%% Problem 5, Part (d.)
% Defining the "truth" relative state vectors over time
rho_truth_RSW = rho_chas_mat_RSW;

% Mean Motion
n = sqrt(mu_earth/sma^3); % Hertz

% Defining RSW components of inital relative pos & vel in RSW frame
rho0_R = rho_0_chas_RSW(1); % km
rho0_S = rho_0_chas_RSW(2); % km
rho0_W = rho_0_chas_RSW(3); % km

zeta_rel_RSW = zeta_0_chas_RSW - cross([0;0;n],rho_0_chas_RSW); % km

zeta0_R = zeta_rel_RSW(1); % km/s
zeta0_S = zeta_rel_RSW(2); % km/s
zeta0_W = zeta_rel_RSW(3); % km/s

for ii = 1:length(tspan)
    
    t = tspan(ii);
    
    rho_CW_R(ii) = (4-3*cos(n*t))*rho0_R + sin(n*t)/n*zeta0_R + 2/n*(1-cos(n*t))*zeta0_S;
    rho_CW_S(ii) = 6*(sin(n*t)-n*t)*rho0_R + rho0_S +2/n*(cos(n*t)-1)*zeta0_R + 1/n*(4*sin(n*t)-3*n*t)*zeta0_S;
    rho_CW_W(ii) = cos(n*t)*rho0_W + sin(n*t)/n*zeta0_W;
    
%     zeta_CW_R(ii) = 3*n*sin(n*t)*rho0_R + cos(n*t)*zeta0_R + 2*sin(n*t)*zeta0_S;
%     zeta_CW_S(ii) = 6*n*(cos(n*t)-1)*rho0_R - 2*sin(n*t)*zeta0_R + (4*cos(n*t)-3)*zeta0_S;
%     zeta_CW_W(ii) = -n*sin(n*t)*rho0_W + cos(n*t)*zeta0_W;
    
end

% Compiling the CW-calculated relative position RSW components
rho_CW_RSW = [rho_CW_R; rho_CW_S; rho_CW_W]; % km

% Error between truth & CW relative position vecs over time in RSW frame
error_rho_RSW = rho_truth_RSW - rho_CW_RSW; % km

% Define error for each component
error_R = error_rho_RSW(1,:); % km
error_S = error_rho_RSW(2,:); % km
error_W = error_rho_RSW(3,:); % km

% Plot the errors on 3 subplots
figure
hold on
subplot(3,1,1)
plot(tspan/T3*3,error_R*1000)
grid on
title('5d) R-Component Error between Truth & CW Relative Position Over Three Orbit Periods')
xlabel('Time (# of Orbit Periods)')
ylabel('R-Component Error (m)')

subplot(3,1,2)
plot(tspan/T3*3,error_S*1000)
grid on
title('5d) S-Component Error between Truth & CW Relative Position Over Three Orbit Periods')
xlabel('Time (# of Orbit Periods)')
ylabel('S-Component Error (m)')

subplot(3,1,3)
plot(tspan/T3*3,error_W*1000)
grid on
title('5d) W-Component Error between Truth & CW Relative Position Over Three Orbit Periods')
xlabel('Time (# of Orbit Periods)')
ylabel('W-Component Error (m)')

hold off

%% Problem 5, Part (e.)
% Magnitude of error at the final time
error_f_mag = norm(error_rho_RSW(:,length(tspan)))*1000; % meters

fprintf('-------------- PROBLEM 5, PART (e) -------------- \n\n')
fprintf('The Magnitude of the Error at the Final Time in meters is:\n')
disp(error_f_mag)

%%
time_to_run_in_seconds = toc