%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ASE 366L: Applied Orbital Mechanics
% HOMEWORK 8 MAIN FILE
tic
close all; format long g; clc; clear all;
% set(0, 'DefaultAxesFontSize', 10, 'defaultlinelinewidth',1)

%% Problem 1
% Givens
mu_earth = 398600.4415; % km^3/s^2
R_earth = 6378.1363; % radius of earth (km)
R_sun = 432288; % km
r_earth_sun_vec = [149597870 0 0]; % Radius vec between Earth & Sun (km)

sma = 20000; % km
ecc = 0;
inc = pi/2; % rad
RAAN = 0; % rad
argp = 0; % rad
tru = 0; % rad

% Vector of Keplarian Orbital Elements
xKOE = [sma, ecc, inc, RAAN, argp, tru];

% Convert initial KOEs into Cartesian position/velocity
Cart0 = kep2cart_Whitbeck(mu_earth,xKOE);

% Calculate Orbit Period, T, & create matching timespan
T = ceil(2 * pi * sqrt(sma^3 / mu_earth)); % sec
tspan = (0:1:T); % 1 orbit period (sec)
[m,n] = size(tspan);

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
r_sat_mat = Y(1:3,:);
v_sat_mat = Y(4:6,:);

% Unit vector from Earth to Sun
r_earth_sun_unit = r_earth_sun_vec ./ norm(r_earth_sun_vec);

gamma = zeros(1,length(tspan));
count = 0;
for ii = 1:length(tspan)
    
    r_sat_vec = r_sat_mat(:,ii);
    r_sat_mag = norm(r_sat_vec);
    test = -sqrt(r_sat_mag^2 - R_earth^2);
    
    % If statement determing if satellite is in shadow or not
    if dot(r_sat_vec, r_earth_sun_unit) < test
        gamma(ii) = 0; % inside the shadow
        count = count + 1; % # of iterations inside the shadow
    else
        gamma(ii) = 1; % outside the shadow
    end    
    
    % If statement checking to see when satellite enters or exits shadow
    if ii == 1
        continue
    elseif (gamma(ii) - gamma(ii-1)) == 1
        exited_shadow = ii; % sec
    elseif (gamma(ii) - gamma(ii-1)) == -1
        entered_shadow = ii; % sec
    else
        continue
    end
    
end

enter_shadow_arg_lat = entered_shadow / T * 360; % degrees
exit_shadow_arg_lat = exited_shadow / T * 360; % degrees
percent_orbit_in_shadow = count / length(tspan) * 100; % percent

fprintf('-------------- PROBLEM 1 -------------- \n\n')

disp('Argument of Latitude (in degrees) when ENTERING the Earth''s Shadow:')
disp(enter_shadow_arg_lat)

disp('Argument of Latitude (in degrees) when EXITING the Earth''s Shadow:')
disp(exit_shadow_arg_lat)

disp('The percent of the orbit spent in the Earth''s Shadow:')
disp(percent_orbit_in_shadow)

%% Problem 2
r_sat_vec_2a = [-7000 -6400 600]; % km
r_sat_vec_2b = [-1000 4695 4316]; % km
r_earth_vec = [0 0 0]; % km

% Part (a)
[gamma_2a,Area_2a] = calc_gamma_shadow(r_sat_vec_2a,r_earth_sun_vec,R_earth,R_sun);

fprintf('-------------- PROBLEM 2, PART (a) -------------- \n\n')
disp('The Fraction of the Sun Visible to the Satellite:')
disp(gamma_2a)
disp('The Approximate Area of Overlap (in km^2):')
disp(Area_2a)

% Part (b)
[gamma_2b,Area_2b] = calc_gamma_shadow(r_sat_vec_2b,r_earth_sun_vec,R_earth,R_sun);

fprintf('-------------- PROBLEM 2, PART (b) -------------- \n\n')
disp('The Fraction of the Sun Visible to the Satellite:')
disp(gamma_2b)
disp('The Approximate Area of Overlap (in km^2):')
disp(Area_2b)

%% Problem 3

% Position & Velocity Vectors of Target Spacecraft in inertial (ijk) Frame
r_targ_vec = [1.638 -0.199 1.130]; % DU
v_targ_vec = [-0.406 -0.100 0.570]; % DU/TU

% Position & Velocity Magnitudes of Target Spacecraft in inertial (ijk) Frame
r_targ_mag = norm(r_targ_vec); % DU
v_targ_mag = norm(v_targ_vec); % DU/TU

% Position of Chaser Spacecraft in RSW frame
rho_chas_vec = [0.1 -0.1 0.2]'; % DU

% Calculating necessary axes
R_unit = [r_targ_vec / r_targ_mag]';
W_unit = [cross(r_targ_vec,v_targ_vec) ./ (r_targ_mag * v_targ_mag)]';
S_unit = cross(W_unit,R_unit);

% Defining Transformation Matrix for RSW -> ijk
T_RSW_ijk = [R_unit S_unit W_unit];

% Calculate relative position of chaser (to target) in ijk frame
r_chas_rel_vec = T_RSW_ijk * rho_chas_vec;

% Calculate absolute position of chaser in ijk frame
r_chas_abs_vec = r_targ_vec' + r_chas_rel_vec;

fprintf('-------------- PROBLEM 3 -------------- \n\n')
disp('The absolute position of the chaser in the i-j-k (inertial) frame (in DU) is:')
disp(r_chas_abs_vec)

%% Problem 4
% Givens
r_vec = [2781 5318 -5629]'; % Pos of s/c (3x1 vec), km
r_earth_sun = [2379260 148079334 -1009936]'; % Pos of Sun wrt Earth (3x1 vec), km
gam = 1; % Ignoring Earth's shadow
C_R = 1.5; % Coef of Reflectivity
A2m = 0.01; % Area-to-mass ratio, m^2/kg

[a_SRP] = accel_SRP(r_vec,r_earth_sun, gam, C_R, A2m);

fprintf('-------------- PROBLEM 4 -------------- \n\n')
disp('The Acceleration of the Spacecraft due to SRP (in km/s^2) is:')
disp(a_SRP)

%% Problem 5
% KOE Givens
SMA = 6800.0; % km
ecc = 0.005;
inc = 71*pi/180; % rad
RAAN = 300*pi/180; % rad
argp = 78*pi/180; % rad
tru = 0*pi/180; % rad
KOE0 = [SMA,ecc,inc,RAAN,argp,tru];

% Other Givens
mu_sun = 1.3271244e11; % km^3/s^2
t_0 = 2458200.5; % Initial time in UTC Julian Days
dUT1 = 0; % UTC to UT1 conversion factor
sec2day = 1/3600/24;% Conversion factor for seconds to days
C_D = 2.0; % Coefficient of Drag
A2m = 0.01; % Area-to-Mass Ratio (m^2/kg)
gam = 1; % Fraction of Sun Visible (ignoring Earth's shadow)
C_R = 1.5; % Coef of Reflectivity
AU2km = 149597870.7; % Conversion factor for AU->km

% Timespan for 24 hours
tspan = (0 : 1 : 24*3600); % sec
[m,n] = size(tspan);

% Convert initial KOEs into Cartesian position/velocity
[Cart0] = kep2cart_Whitbeck( mu_earth, KOE0 );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ODE45 Setup %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ODE options
odeoptions = odeset('RelTol',1e-25,'AbsTol',1e-30);

%% Problem 5, Part (a)
% Anonymous function outputting 6x1 vector w/ 3D velocity & acceleration
% w/ input of time & 6x1 3D position/velocity vector.
vel_acc_5a = @(t,Cart) [Cart(4:6);...
    -mu_earth .* Cart(1:3) ./ (norm(Cart(1:3))^3)...
    + accel_J2_J3(Cart(1:3))...
    + accel_TB(mu_sun,Cart(1:3),t,t_0,dUT1)...
    + accel_drag(Cart,C_D,A2m)...
    + accel_SRP(Cart(1:3),AU2km*sun(t*sec2day+t_0,dUT1), gam, C_R, A2m)];

% Calling ode45 function to output time vector T & position/velocity matrix
% at each time interval Y.
[~,Y_5a] = ode45(vel_acc_5a, tspan, Cart0, odeoptions);

% Get Y matrix in form where columns contain pos/vel & each new column is
% another time t
Y_5a = Y_5a';
r_t_mat_5a = Y_5a(1:3,:);
v_t_mat_5a = Y_5a(4:6,:);

% Final position & velocity vectors
r_f_vec = r_t_mat_5a(:,n);
v_f_vec = v_t_mat_5a(:,n);

fprintf('-------------- PROBLEM 5, PART (a) -------------- \n\n')
disp('The Final Position Vector (in km) is:')
disp(r_f_vec)
disp('The Final Velocity Vector (in km/s) is:')
disp(v_f_vec)

%% Problem 5, Part (b)
% Anonymous function outputting 6x1 vector w/ 3D velocity & acceleration
% w/ input of time & 6x1 3D position/velocity vector.
vel_acc_5b = @(t,Cart) [Cart(4:6);...
    -mu_earth .* Cart(1:3) ./ (norm(Cart(1:3))^3)];

% Calling ode45 function to output time vector T & position/velocity matrix
% at each time interval Y.
[~,Y_5b] = ode45(vel_acc_5b, tspan, Cart0, odeoptions);

% Get Y matrix in form where columns contain pos/vel & each new column is
% another time t
Y_5b = Y_5b';
r_t_mat_5b = Y_5b(1:3,:);
v_t_mat_5b = Y_5b(4:6,:);

% Calc difference between positions for parts (b) & (a) over time
dr_dt = r_t_mat_5b - r_t_mat_5a;
dr_i = dr_dt(1,:); % i-component of position over time
dr_j = dr_dt(2,:); % j-component of position over time
dr_k = dr_dt(3,:); % k-component of position over time

% Plot the differences
figure
hold on

subplot(3,1,1)
plot(tspan/3600,dr_i)
title('Difference in r_i between 2-Body-only & Additional Forces over Time')
xlabel('Time (hr)')
ylabel('Difference in r_i (km)')
grid on

subplot(3,1,2)
plot(tspan/3600,dr_j)
title('Difference in r_j between 2-Body-only & Additional Forces over Time')
xlabel('Time (hr)')
ylabel('Difference in r_j (km)')
grid on

subplot(3,1,3)
plot(tspan/3600,dr_k)
title('Difference in r_k between 2-Body-only & Additional Forces over Time')
xlabel('Time (hr)')
ylabel('Difference in r_k (km)')
grid on

hold off

%% Problem 5, Part(c)
% Calc the magnitude of the final position error
dr_final_mag = norm(dr_dt(:,n));

fprintf('-------------- PROBLEM 5, PART (c) -------------- \n\n')
disp('The Magnitude of the Position Error (in km) is:')
disp(dr_final_mag)


%%
time_to_run_in_seconds = toc