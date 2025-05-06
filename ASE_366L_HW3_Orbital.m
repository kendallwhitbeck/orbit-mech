%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ASE 366L: Applied Orbital Mechanics
% HOMEWORK 3 MAIN FILE

close all; format long g; clc; %clear all;
set(0, 'DefaultAxesFontSize', 10, 'defaultlinelinewidth',1)

%% Prob 2
disp('Problem 2 %%%%%%%%%%%%%%%%%')
% Givens:
year = 2018;
month = 2;
day = 9;
hour = 15;
min = 0;
sec = .184798;

JD_UT1 = 367*year - floor(7*(year + floor((month+9)/12))/4) +floor(275*month/9)...
    + day + 1721013.5 + 1/24*(hour + 1/60*(min + sec/60)); % Days

MJD_UT1 = JD_UT1 - 2400000.5; % Days

T_UT1 = (MJD_UT1 - 2451545 + 2400000.5)/36525; % centuries!

GMST_sec = 67310.54841 + (876600*3600 + 8640184.812866)*T_UT1 + ...
    .093104*T_UT1^2 - 6.2e-6 * T_UT1^3; % sec

theta_GMST_deg_total = GMST_sec/3600 * 15; % degrees total

theta_GMST_deg_360_2 = mod(theta_GMST_deg_total,360); % Degrees btwn 0 & 360

%% Prob 3
disp('Problem 3 %%%%%%%%%%%%%%%%%')

[JD,MJD,T,GMST_sec,theta_GMST_deg_total,theta_GMST_deg_360] = calc_GMST_ang(2016,10,24,6,30,30);

%% Prob 4
disp('Problem 4 %%%%%%%%%%%%%%%%%')
% clear all; clc
% Givens:
mu = 398600.4415; % km^3/s^2
dAT = 37; % seconds
t_p = 0; % seconds

r_vec_A = [-1140.697 6126.252 -926.199]; % km
v_vec_A = [-5.706065 -1.934816 -5.7701134]; % km/s
UTC_A = 2458165.4375; % JD format (days)

r_vec_B = [-1350.501 + 6048.007 -1138.605]; % km
v_vec_B = [-5.632593 -2.293872 -5.709216]; % km/s
TT_B = 2458165.4375;% JD format (days)

% Convert TT_B to UTC
TAI_B = TT_B - 32.184; % seconds
UTC_B = TAI_B - dAT; % seconds

% Diffence from UTC_A to UTC_B
dUTC = UTC_B - UTC_A;

% Convert Cart2Kep at TT_A
[ KOE_TT_A ] = cart2kep_Whitbeck(mu, [r_vec_A, v_vec_A]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Need satellite A pos/vel at dTT seconds ago from TT_A
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculates Mean Anomaly for elliptical orbit
Me = sqrt(mu/KOE_TT_A(1)^3)*(dUTC-t_p);

% Ensure that MEan Anomaly remains under 2*pi radians
while Me > 2*pi
    Me = Me - 2*pi;
end

% Call function converting mean anomaly to true anomaly
[tru_TT_B] = Me2tru(KOE_TT_A(2),Me); % radians

% Assign new value of true anomaly (tru_TT_B)
KOE_TT_B = KOE_TT_A;
KOE_TT_B(6) = tru_TT_B;

% Kep2Cart at time TT_B to get r and v
[ xC ] = kep2cart_Whitbeck( mu, KOE_TT_B )
fprintf('The two operators are not tracking the same satellite.\n\n')

%% Prob 5
disp('Problem 5 %%%%%%%%%%%%%%%%%')

[JD_UT1,theta_ERA,W,R,PN,T_ITRF2GCRF] = ...
    ITRF2GCRF(2018,02,9,15,0,0,.00575134,.30394248,.184798,0,0,358.4586,-7.138811,.004737961)

%% Prob 6
disp('Problem 6 %%%%%%%%%%%%%%%%%')

r_ITRF_vec = [-742.845 -5463.244 3196.066]'; % km

disp(' -Part a:')
r_GCRF_vec_a = calc_R3(-theta_GMST_deg_360_2*pi/180) * r_ITRF_vec % km

disp(' -Part b:')
r_GCRF_vec_b = T_ITRF2GCRF * r_ITRF_vec % km

disp(' -Part c:')
error = r_GCRF_vec_a - r_GCRF_vec_b % km
