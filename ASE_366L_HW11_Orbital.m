%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ASE 366L: Applied Orbital Mechanics
% HOMEWORK 11 MAIN FILE
tic
close all; format long g; clc; clear all;
% set(0, 'DefaultAxesFontSize', 10, 'defaultlinelinewidth',1)

%% Problem 1
% Givens
R_Earth = 6378.1363; % km
w_Earth = 7.2927e-5; % Angular Velocity of Earth in rad/s
mu_Earth = 398600.4415; % km^3/s^2
% time coordinates
t1 = 0; % sec
t2 = 20; % sec
% Range Observation (p) coords
p1 = 3590.0; % km
p2 = 3610.0; % km
% Azimuth Observation (az) coords
az1 = 15.00*pi/180; % rad
az2 = 14.19*pi/180; % rad
% Elevation Observation (el) coords
el1 = 76.09*pi/180; % rad
el2 = 74.15*pi/180; % rad

% Anonymous Function for Relative Observation Position Vector in SEZ frame
p_SEZ_vec = @(p,az,el) [-p*cos(el)*cos(az), p*cos(el)*sin(az), p*sin(el)]';

% Observation Position Vectors in SEZ frame
p1_SEZ_vec = p_SEZ_vec(p1, az1, el1); % km
p2_SEZ_vec = p_SEZ_vec(p2, az2, el2); % km

fprintf('-------------- PROBLEM 1, PART (a) -------------- \n\n')
fprintf('Observation Position Vector in SEZ frame at t=0s in km is:\n')
disp(p1_SEZ_vec)
fprintf('Observation Position Vector in SEZ frame at t=20s in km is:\n')
disp(p2_SEZ_vec)

%% Problem 1, Part (b)
% Latitude & Longitude are zero bc R_vec = R_earth in i-direction only
lat = 30*pi/180; % latitude (rad)
long = -98*pi/180; % longitude (rad)

% Transformation matrix from SEZ to ITRF
T_SEZ_ITRF = [calc_R2(pi/2-lat) * calc_R3(long)]';
T_SEZ_ITRF(abs(T_SEZ_ITRF)<1e-15) = 0;

fprintf('-------------- PROBLEM 1, PART (b) -------------- \n\n')
fprintf('The Rotation Matrix from SEZ frame to ITRF frame is:\n')
disp(T_SEZ_ITRF)

%% Problem 1, Part (c)
% Transformation matrix from ITRF to GCRF for each time
T1_ITRF_GCRF = calc_R3(-w_Earth*t1);
T2_ITRF_GCRF = calc_R3(-w_Earth*t2);

% Observer's position vector for location on Earth
R_ITRF_vec = [R_Earth*cos(lat)*cos(long)
              R_Earth*cos(lat)*sin(long)
              R_Earth*sin(lat)]; % km

% Position vector of object in GCRF frame at each time
r1_GCRF_vec = T1_ITRF_GCRF * (T_SEZ_ITRF * p1_SEZ_vec + R_ITRF_vec);
r2_GCRF_vec = T2_ITRF_GCRF * (T_SEZ_ITRF * p2_SEZ_vec + R_ITRF_vec);

fprintf('-------------- PROBLEM 1, PART (c) -------------- \n\n')
fprintf('The Object''s Position Vector in GCRF frame at t=0s in km is:\n')
disp(r1_GCRF_vec)
fprintf('The Object''s Position Vector in GCRF frame at t=20s in km is:\n')
disp(r2_GCRF_vec)

%% Problem 1, Part (d)
% Change in time from t1 to t2
dt = t2 - t1; % sec

% Assume short way transfer s.t. tm = +1, not -1
tm = 1;

[v1_vec,v2_vec] = Gauss_Lambert(r1_GCRF_vec,r2_GCRF_vec,dt,tm,mu_Earth);

fprintf('-------------- PROBLEM 1, PART (d) -------------- \n\n')
fprintf('The Object''s Velocity Vector in GCRF frame at t=0s in km/s is:\n')
disp(v1_vec)
fprintf('The Object''s Velocity Vector in GCRF frame at t=20s in km/s is:\n')
disp(v2_vec)

%% Problem 2, Part (a)
% Givens
r1_vec2 = [ 2 1]'; % DU
r2_vec2 = [-2 1]'; % DU
mu = 1; % DU^3/TU^2

% Magnitude of position vectors
r1 = norm(r1_vec2); % DU
r2 = norm(r2_vec2); % DU

% Cosine of change in true anomaly
cos_dv = dot(r1_vec2,r2_vec2) / (r1*r2);

% Chord line from r1_vec to r2_vec
c = sqrt(r1^2 + r2^2 - 2*r1*r2*cos_dv); % DU

% Semiperimeter
s = (r1+r2+c)/2; % DU

% Minimum semimajor axis for the orbit that connects r1_vec & r2_vec
a_min = s/2; % DU

fprintf('-------------- PROBLEM 2, PART (a) -------------- \n\n')
fprintf('The Minimum SMA for the Orbit connecting r1 & r2 in DU is:\n')
disp(a_min)

%% Problem 2, Part (b)
sma = 2.5; % DU

R1 = 2*sma - r1;
R2 = 2*sma - r2;

dy1 = sqrt(R1^2 - r1_vec2(1)^2) + r1_vec2(2);
dy2 = sqrt(R2^2 - r1_vec2(1)^2) - r1_vec2(2);

ecc1 = dy1 / (2*sma);
ecc2 = dy2 / (2*sma);

fprintf('-------------- PROBLEM 2, PART (b) -------------- \n\n')
fprintf('The Eccentricity for the First Orbit connecting r1 & r2 is:\n')
disp(ecc1)
fprintf('The Eccentricity for the Second Orbit connecting r1 & r2 is:\n')
disp(ecc2)

%% Problem 3
mu_Sun = 1; % AU^3/TU^2

% Earth's position vector at time of departure in heliocentric frame
r_Earth = [1 0 0]'; % AU

% Mar's position vector at time of arrival in heliocentric frame
r_Mars = [0 1.524 0]'; % AU

[v1_vec,t_abs,a_min,t_0f_min] = Lambert_E_min(r_Earth,r_Mars,mu_Sun);

fprintf('-------------- PROBLEM 3, PART (a) -------------- \n\n')
fprintf('The Absolute Minimum Time from Earth to Mars in TU is:\n')
disp(t_abs)
fprintf('-------------- PROBLEM 3, PART (b) -------------- \n\n')
fprintf('The Minimum SMA from Earth to Mars in DU is:\n')
disp(a_min)
fprintf('-------------- PROBLEM 3, PART (c) -------------- \n\n')
fprintf('The velocity required (heliocentric frame) for transfer at departure in DU/TU is:\n')
disp(v1_vec)
fprintf('-------------- PROBLEM 3, PART (d) -------------- \n\n')
fprintf('The Short Way Transfer Time using Minimum Energy in TU is:\n')
disp(t_0f_min)


%% Problem 4
mu_Sun = 1; % AU^3/TU^2

% Earth's position vector at time of departure in heliocentric frame
r_Earth = [1 0 0]'; % AU

% Mar's KOEs
tru_vec = (5:5:170)*pi/180; % rads
for ii = 1:length(tru_vec)
    xKOE_mars = [1.524, 0, 1.85*pi/180, 0, 0, tru_vec(ii)]';
    r_Mars4 = kep2cart_Whitbeck( mu_Sun, xKOE_mars );
    [v1_vec4(ii,:),t_abs4(ii),a_min4(ii),t_0f_min4(ii)] = Lambert_E_min(r_Earth,r_Mars4(1:3),mu_Sun);
end

%% Problem 4, Part (a)
figure
hold on
plot(tru_vec*180/pi,t_abs4)
grid on
title('4a. Absolute Minimum Time for Transfer as a Function of Mar''s True Anomaly in Heliocentric Frame')
xlabel('Mar''s True Anomaly in Heliocentric Frame (degrees)')
ylabel('Absolute Minimum Time (TU)')
hold off

%% Problem 4, Part (b)
sma_Earth = 2*r_Earth(1)/2; % AU
KOE_Earth = [sma_Earth 0 0 0 0 0]';
Cart_Earth = kep2cart_Whitbeck(mu_Sun,KOE_Earth);
vel_Earth = Cart_Earth(4:6);

fprintf('-------------- PROBLEM 4, PART (b) -------------- \n\n')
fprintf('The Velocity Vector of Earth (Heliocentric Frame) in AU/TU is:\n')
disp(vel_Earth)

%% Problem 4, Part (c)
for ii = 1:length(tru_vec)
    v_inf_vec = norm(v1_vec4(ii,:) - vel_Earth');
    C3(ii) = v_inf_vec^2;
end

figure
hold on
plot(tru_vec*180/pi,C3)
grid on
title('4c. Departure C_3 for Transfer as a Function of Mar''s True Anomaly in Heliocentric Frame')
xlabel('Mar''s True Anomaly in Heliocentric Frame (degrees)')
ylabel('Departure C_3 (AU^2/TU^2)')
hold off


%%
toc