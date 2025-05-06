clc
clear all
format shortg

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 2 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Given:
mu = 398600.4415; % Standard Gravitational Parameter of Earth in km^3/s^2
r_vec = [-17130, 39, -784]; % km
v_vec = [-.251, -2.827, 3.754]; % km/s
t = 0; % sec

[h_vec r epsilon e_vec e theta E M_e a r_p r_a p t_p] = ...
    orbit_parameter_calculation(r_vec,v_vec,mu,t);

%% Problem 7 %

% time vector for plotting
t_vec = 0:500:40000; % seconds

% Calculates true anomaly, theta, as a function of time
theta_vec = mu^2 / norm(h_vec).^3 .* (t_vec - t_p);
for i = 1:length(theta_vec)
    theta_vec(i) = mod(theta_vec(i),2*pi);
    
    if t_vec(i) == 15000
        theta_15000s = theta_vec(i)
    end
end

% Calculates postion vector, r_vec, as a function of time
r_vec = p ./ (1 + e .* cos(theta_vec));

% Calculates orbital velocity vector, v_vec, as a function of time
v_vec = sqrt(mu * ((2 ./ r_vec) - (1 / a)));

% Calculates total orbital engergy of initial conditions
epsilon_vec = (v_vec.^2 ./ 2) - (mu ./ r_vec);
epsilon_vec = epsilon_vec - epsilon;

%% Plots theta, r and v as functions time
figure
set(0, 'DefaultAxesFontSize', 12, 'DefaultLineLineWidth', 2)

hold on

subplot(3,1,1)
plot(t_vec,theta_vec)
title('True Anomaly over Time for the Given Orbit of a Spacecraft Revolving around Earth')
xlabel('Time (s)')
ylabel('True Anomaly (radians)')

subplot(3,1,2)
plot(t_vec,r_vec)
title('Position over Time for the Given Orbit of a Spacecraft Revolving around Earth')
xlabel('Time (s)')
ylabel('Position (km)')

subplot(3,1,3)
plot(t_vec,v_vec)
title('Velocity over Time for the Given Orbit of a Spacecraft Revolving around Earth')
xlabel('Time (s)')
ylabel('Velocity (km/s)')

figure
plot(t_vec,epsilon_vec)
title('Change in Specific Energy over Time for the Given Orbit of a Spacecraft Revolving around Earth')
xlabel('Time (s)')
ylabel('Change in Specific Energy (km^2/s^2)')

hold off





