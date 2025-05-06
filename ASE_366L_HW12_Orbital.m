%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ASE 366L: Applied Orbital Mechanics
% HOMEWORK 12 MAIN FILE
tic
close all; format long g; clc; %clear all;
set(0, 'DefaultAxesFontSize', 16, 'defaultlinelinewidth',3)

%% Problem 1, Part (a)
% Standard Gravitational Parameter of the Sun
mu_Sun = 1.32712440018e11; % km^3/s^2
% Earth's orbital elements in heliocentric frame (km & rads)
KOE_Earth = [149599802, .0167, 0, 174.87*pi/180, -71.935*pi/180, 0]';
Cart_Earth = kep2cart_Whitbeck(mu_Sun,KOE_Earth)';
pos_Earth = Cart_Earth(1:3)';
vel_Earth = Cart_Earth(4:6)';
% Range for True Anomaly of Mars
tru_vec = (1:1:359)*pi/180; % rads
tm = 1;
for ii = 1:length(tru_vec)
    % Mar's orbital elements in heliocentric frame (km & rads)
    KOE_Mars = [227941896, .093, 1.84*pi/180, 49.558*pi/180, 286.5*pi/180, tru_vec(ii)]';
    
    % Mar's Cartesian State at each iteration (heliocentric frame)
    Cart_Mars = kep2cart_Whitbeck( mu_Sun, KOE_Mars ); % km & km/s
    pos_Mars = Cart_Mars(1:3); % position in km
    vel_Mars = Cart_Mars(4:6); % Velocity in km/s
    
    % Finds departure vel & arrival vel needed for transfer orbit to Mars
    [v1_vec,v2_vec,~,~,~] = Lambert_E_min(pos_Earth,pos_Mars,mu_Sun);
    
    % Excess Hyperbolic Velocity Needed to Leave Earth
    v_inf_Earth = norm(v1_vec - vel_Earth);
    % Burn Energy Needed to Escape Earth on desire Trajectory
    C3(ii) = v_inf_Earth^2;
    
    % Hyperbolic Excess Velocity when Arriving at Mars
    v_inf_Mars(ii) = norm(v2_vec - vel_Mars);
end
figure
hold on
subplot(2,1,1)
plot(tru_vec*180/pi,C3)
grid on
title('1a. Departure C_3 for Transfer vs. \nu_{Mars} in Heliocentric Frame')
xlabel('\nu_{Mars} in Heliocentric Frame (degrees)')
ylabel('Departure C_3 (km^2/s^2)')

subplot(2,1,2)
plot(tru_vec*180/pi,v_inf_Mars)
grid on
title('1a. Hyperbolic Excess Velocity when Arriving at Mars vs. \nu_{Mars} in Heliocentric Frame')
xlabel('\nu_{Mars} in Heliocentric Frame (degrees)')
ylabel('v_{\infty,Mars} (km/s)')
hold off

%% Problem 1, Part (b.)
mindex = find(C3 == min(C3(:))); % index # containing minimum C3 value
tru_minC3 = tru_vec(mindex)*180/pi; % deg

fprintf('-------------- PROBLEM 1, PART (b) -------------- \n\n')
fprintf('The True Anomaly of Mars w/ Minimum C3 value (in deg) is:\n')
disp(tru_minC3)

%% Problem 1, Part (c.)
mindex = find(v_inf_Mars == min(v_inf_Mars(:))); % index # containing minimum C3 value
tru_minVinf = tru_vec(mindex)*180/pi; % deg

fprintf('-------------- PROBLEM 1, PART (c) -------------- \n\n')
fprintf('The True Anomaly of Mars w/ Minimum Hyperbolic \nExcess Velocity at Mars Arrival (in deg) is:\n')
disp(tru_minVinf)

%% Problem 2
clearvars -except mu_Sun

% Earth's orbital elements in heliocentric frame (km & rads)
KOE_Earth = [149599802, .0167, 0, 174.87*pi/180, -71.935*pi/180, 0]';
% Mar's orbital elements in heliocentric frame (km & rads)
KOE_Mars = [227941896, .093, 1.84*pi/180, 49.558*pi/180, 286.5*pi/180, 0]';
% True Anomaly Range for both
tru_vec2 = (0:1:359)*pi/180; % rad

for E = 1:length(tru_vec2)
    for M = 1:length(tru_vec2)

        % Update True Anomaly for each
        KOE_Mars(6) = tru_vec2(M);
        KOE_Earth(6) = tru_vec2(E);

        % Mar's Cartesian State at each iteration (heliocentric frame)
        Cart_Mars = kep2cart_Whitbeck( mu_Sun, KOE_Mars ); % km & km/s
        pos_Mars = Cart_Mars(1:3); % position in km
        vel_Mars = Cart_Mars(4:6); % Velocity in km/s
        
        % Earth's Cartesian State at each iteration (heliocentric frame)
        Cart_Earth = kep2cart_Whitbeck(mu_Sun,KOE_Earth)';
        pos_Earth = Cart_Earth(1:3)';
        vel_Earth = Cart_Earth(4:6)';

        % Finds departure vel & arrival vel needed for transfer orbit to Mars
        [v1_vec,v2_vec,~,~,t_0f_min(M,E)] = Lambert_E_min(pos_Earth,pos_Mars,mu_Sun);
        
        % Excess Hyperbolic Velocity Needed to Leave Earth
        v_inf_Earth = norm(v1_vec - vel_Earth);
        % Burn Energy Needed to Escape Earth on desire Trajectory
        C3_2(M,E) = v_inf_Earth^2;
        
        % Hyperbolic Excess Velocity when Arriving at Mars
        v_inf_Mars_2(M,E) = norm(v2_vec - vel_Mars);
        
    end
end
% Converting t_0f from seconds to days
t_0f_min = t_0f_min /3600/24; % days

%% Problem 2, Part (a.)
figure
hold on
contour(tru_vec2*180/pi,tru_vec2*180/pi,C3_2)
grid on
title('2a. Departure C_3 for Earth-Mars Transfer Opportunities')
xlabel('\nu_{Earth} in Heliocentric Frame (degrees)')
ylabel('\nu_{Mars} in Heliocentric Frame (degrees)')
c = colorbar;
c.Label.String = 'Departure C_3 (km^2/s^2)';
hold off

%% Problem 2, Part (b.)
figure
hold on
contour(tru_vec2*180/pi,tru_vec2*180/pi,v_inf_Mars_2)
grid on
title('2b. Hyperbolic Excess Velocity at Arrival (v_{\infty,Mars}) for Earth-Mars Transfer Opportunities')
xlabel('\nu_{Earth} in Heliocentric Frame (degrees)')
ylabel('\nu_{Mars} in Heliocentric Frame (degrees)')
c = colorbar;
c.Label.String = 'v_{\infty,Mars} (km/s)';
hold off

%% Problem 2, Part (c.)
figure
hold on
contour(tru_vec2*180/pi,tru_vec2*180/pi,t_0f_min)
grid on
title('2c. Minimum Energy Time (t_{of}) for Earth-Mars Transfer Opportunities')
xlabel('\nu_{Earth} in Heliocentric Frame (degrees)')
ylabel('\nu_{Mars} in Heliocentric Frame (degrees)')
c = colorbar;
c.Label.String = 't_{of} (days)';
hold off

%% Problem 2, Part (d.)
% Min value of C3
min_C3_2 = min(C3_2(:));

% Index location of min C3 value
[index_Mars_C3_min,index_Earth_C3_min] = find(C3_2==min_C3_2);

% Converting index values to degrees
v_Earth_C3_min = index_Earth_C3_min - 1; % deg
v_Mars_C3_min = index_Mars_C3_min - 1; % deg

fprintf('-------------- PROBLEM 2, PART (d) -------------- \n\n')
fprintf('The True Anomaly Combo of Earth & Mars w/ Minimum C3 value is:\n\n')
fprintf('Minimum C3 value at Departure (in km^2/s^2) is:\n')
disp(min_C3_2)
fprintf('True Anomaly of Earth (in deg)\n')
disp(v_Earth_C3_min)
fprintf('True Anomaly of Mars (in deg)\n')
disp(v_Mars_C3_min)

%% Problem 2, Part (e.)
% Min value of C3
min_v_inf_Mars_2 = min(v_inf_Mars_2(:));

% Index location of min C3 value
[index_Mars_v_inf_Mars_2_min,index_Earth_v_inf_Mars_2_min] = find(v_inf_Mars_2==min_v_inf_Mars_2);

% Converting index values to degrees
v_Earth_v_inf_Mars_2_min = index_Earth_v_inf_Mars_2_min - 1; % deg
v_Mars_v_inf_Mars_2_min = index_Mars_v_inf_Mars_2_min - 1; % deg

fprintf('-------------- PROBLEM 2, PART (e) -------------- \n\n')
fprintf('The True Anomaly Combo of Earth & Mars w/ Minimum \nHyperbolic Excess Velocity at Arrival is:\n\n')
fprintf('Minimum Hyperbolic Excess Velocity at Arrival (in km/s)\n')
disp(min_v_inf_Mars_2)
fprintf('True Anomaly of Earth (in deg)\n')
disp(v_Earth_v_inf_Mars_2_min)
fprintf('True Anomaly of Mars (in deg)\n')
disp(v_Mars_v_inf_Mars_2_min)

%% Problem 2, Part (f.)
% Approximate change in velocity needed for departure
dv1 = C3_2.^.5;
% Approximate change in velocity needed at arrival
dv2 = v_inf_Mars_2;
% calc total change in velocity
dv_tot = abs(dv1) + abs(dv2);
% Min index locations
[index_M,index_E] = find(dv_tot == min(dv_tot(:)));
% Convert index to true anomaly
tru_E = index_E - 1;
tru_M = index_M - 1;


%%
toc