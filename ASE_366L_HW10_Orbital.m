%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ASE 366L: Applied Orbital Mechanics
% HOMEWORK 10 MAIN FILE
tic
close all; format long g; clc; clear all;
% set(0, 'DefaultAxesFontSize', 10, 'defaultlinelinewidth',1)

%% Problem 1
R_Earth = 6378.1363; % km
w_Earth = 7.2927e-5; % Angular Velocity of Earth in rad/s
sma = 42000; % km
ecc = 0;

% time coordinates
t1 = 0; % sec
t2 = 10*60; % sec

% Right Ascension (RA) coords
ra1 = 18.661*pi/180; % rad
ra2 = 21.547*pi/180; % rad

% Declination (Dec) coords
dec1 = -8.6916*pi/180; % rad
dec2 = -7.9973*pi/180; % rad

% Observer's position vector in ITRF frame
R_GCRF_vec = [R_Earth 0 0]'; % km
R_GCRF_unit = [1 0 0]';
R_GCRF = R_Earth; % km

% Position unit vector of sat at both times
p1_unit = [cos(dec1)*cos(ra1), cos(dec1)*sin(ra1), sin(dec1)]';
p2_unit = [cos(dec2)*cos(ra2), cos(dec2)*sin(ra2), sin(dec2)]';

% Angle btwn R and p
theta1 = pi - arccos(dot(R_GCRF_unit,p1_unit));
theta2 = pi - arccos(dot(R_GCRF_unit,p2_unit));

% Quadratically solve for range magnitude at both times, p1 and p2
coef_vec1 = [1, -2*R_GCRF*cos(theta1), R_GCRF^2-sma^2]';
coef_vec2 = [1, -2*R_GCRF*cos(theta2), R_GCRF^2-sma^2]';
p1 = max(roots(coef_vec1));
p2 = max(roots(coef_vec2));

% Calc absolute pos. unit vector from center of Earth @ t1 & t2
r1_unit = R_GCRF/sma.*R_GCRF_unit + p1/sma.*p1_unit;
r2_unit = R_GCRF/sma.*R_GCRF_unit + p2/sma.*p2_unit;

% Angular Momentum unit vector
h_unit = cross(r1_unit,r2_unit)./norm(cross(r1_unit,r2_unit));

% Inclination
inc = arccos(h_unit(3));

% Right Ascension of the Ascending Node
RAAN = atan2(h_unit(1), -h_unit(2)); % rad

n_unit = [cos(RAAN), sin(RAAN), 0]';

% True Anomaly
u1 = arccos(dot(n_unit,r1_unit))*180/pi; % deg
if r1_unit(3) <0
    u1 = 360 - u1;
end

% Converting to degrees
RAAN = RAAN*180/pi;
inc = inc*180/pi;

fprintf('-------------- PROBLEM 1 -------------- \n\n')
fprintf('The Inclination in degrees is:\n')
disp(inc)
fprintf('The Right Ascension of Ascending Nodes in degrees is:\n')
disp(RAAN)
fprintf('The True Anomaly at t=0s in degrees is:\n')
disp(u1)



%% Problem 2
% time coordinates
t1 = 100*(1-1); % seconds
t2 = 100*(2-1); % seconds
t3 = 100*(3-1); % seconds

%% Problem 2, Part (a.)
% Initial Position Vectors in GCRF Frame
r1_vec = [34159 23915 -417]'; % km
r2_vec = [33979 24167 -422]'; % km
r3_vec = [33796 24418 -426]'; % km

fprintf('-------------- PROBLEM 2, PART (a) -------------- \n\n')

% Using IOD to calculate velocity vector v_2 at t_2
[v2_vec] = IOD(r1_vec,r2_vec,r3_vec,t1,t2,t3);

fprintf('\nThe velocity vector at t=100s is:\n')
disp(v2_vec)

%% Problem 2, Part (b.)
% Initial Position Vectors in GCRF Frame
r1_vec = [-1970 4033 -5413]'; % km
r2_vec = [-2111 3421 -5800]'; % km
r3_vec = [-2228 2770 -6122]'; % km

fprintf('-------------- PROBLEM 2, PART (b) -------------- \n\n')

% Using IOD to calculate velocity vector v_2 at t_2
[v2_vec] = IOD(r1_vec,r2_vec,r3_vec,t1,t2,t3);

fprintf('\nThe velocity vector at t=100s is:\n')
disp(v2_vec)

%% Problem 2, Part (c.)
% Initial Position Vectors in GCRF Frame
r1_vec = [20502 15040 -5516]'; % km
r2_vec = [20418 15263 -5206]'; % km
r3_vec = [20330 15483 -4895]'; % km

fprintf('-------------- PROBLEM 2, PART (c) -------------- \n\n')

% Using IOD to calculate velocity vector v_2 at t_2
[v2_vec] = IOD(r1_vec,r2_vec,r3_vec,t1,t2,t3);

fprintf('\nThe velocity vector at t=100s is:\n')
disp(v2_vec)

%% Problem 3
% Time Observations
t1 = 0; % sec
t2 = 30; % sec
t3 = 60; % sec

% Range Observations
p1 = 651.343; % km
p2 = 865.398; % km
p3 = 1088.356; % km

% Azimuth Observations
az1 = 45.000*pi/180; % rad
az2 = 45.715*pi/180; % rad
az3 = 46.102*pi/180; % rad

% Elevation Observations
el1 = 26.411*pi/180; % rad
el2 = 17.782*pi/180; % rad
el3 = 12.210*pi/180; % rad

% Observer's position vector in ITRF frame
R_ITRF_vec = [R_Earth 0 0]'; % km

%% Problem 3, Part (a)
% Anonymous Function for Observation Position Vector in SEZ frame
p_SEZ_vec = @(p,az,el) [-p*cos(el)*cos(az), p*cos(el)*sin(az), p*sin(el)]';

% Observation Position Vectors in SEZ frame
p1_SEZ_vec = p_SEZ_vec(p1, az1, el1);
p2_SEZ_vec = p_SEZ_vec(p2, az2, el2);
p3_SEZ_vec = p_SEZ_vec(p3, az3, el3);

fprintf('-------------- PROBLEM 3, PART (a) -------------- \n\n')
fprintf('Observation Position Vector in SEZ frame at t=0s in km is:\n')
disp(p1_SEZ_vec)
fprintf('Observation Position Vector in SEZ frame at t=30s in km is:\n')
disp(p2_SEZ_vec)
fprintf('Observation Position Vector in SEZ frame at t=60s in km is:\n')
disp(p3_SEZ_vec)

%% Problem 3, Part (b)
% Latitude & Longitude are zero bc R_vec = R_earth in i-direction only
lat_gd = 0; % Geodetic latitude (rad)
long = 0; % longitude (rad)

% Transformation matrix from SEZ to ITRF
T_SEZ_ITRF = [calc_R2(pi/2-lat_gd) * calc_R3(long)]';
T_SEZ_ITRF(abs(T_SEZ_ITRF)<1e-15) = 0;

fprintf('-------------- PROBLEM 3, PART (b) -------------- \n\n')
fprintf('The Rotation Matrix from SEZ frame to ITRF frame is:\n')
disp(T_SEZ_ITRF)

%% Problem 3, Part (c)
% Transformation matrix from ITRF to GCRF for each time
T1_ITRF_GCRF = calc_R3(-w_Earth*t1);
T2_ITRF_GCRF = calc_R3(-w_Earth*t2);
T3_ITRF_GCRF = calc_R3(-w_Earth*t3);

% Position vector of object in GCRF frame at each time
r1_GCRF_vec = T1_ITRF_GCRF * (T_SEZ_ITRF * p1_SEZ_vec + R_ITRF_vec);
r2_GCRF_vec = T2_ITRF_GCRF * (T_SEZ_ITRF * p2_SEZ_vec + R_ITRF_vec);
r3_GCRF_vec = T3_ITRF_GCRF * (T_SEZ_ITRF * p3_SEZ_vec + R_ITRF_vec);

fprintf('-------------- PROBLEM 3, PART (c) -------------- \n\n')
fprintf('The Object''s Position Vector in GCRF frame at t=0s in km is:\n')
disp(r1_GCRF_vec)
fprintf('The Object''s Position Vector in GCRF frame at t=30s in km is:\n')
disp(r2_GCRF_vec)
fprintf('The Object''s Position Vector in GCRF frame at t=60s in km is:\n')
disp(r3_GCRF_vec)


%% Problem 3, Part (d,e)

fprintf('-------------- PROBLEM 3, PART (d) -------------- \n\n')
% Using IOD to calculate velocity vector v_2 at t_2
[v2_vec] = IOD(r1_GCRF_vec,r2_GCRF_vec,r3_GCRF_vec,t1,t2,t3);

fprintf('\n-------------- PROBLEM 3, PART (e) -------------- \n')
fprintf('\nThe velocity vector v_2 at t_2 is:\n')
disp(v2_vec)

%%
time_to_run_in_seconds = toc