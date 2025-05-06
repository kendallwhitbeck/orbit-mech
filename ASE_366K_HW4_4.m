clc
format shortg

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 4 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Given:
mu = 398600.4415; % Standard Gravitational Parameter of Earth in km^3/s^2
r_vec = [0, 39991, -52979]; % km
v_vec = [3.5511, 0, 0]; % km/s
t = 0; % sec

[h_vec r epsilon e_vec e theta E M_e a r_p r_a p t_p F M_h] = ...
    orbit_parameter_calculation(r_vec,v_vec,mu,t);

%%
% time vector for plotting
t_vec = -20000:500:20000; % seconds

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
epsilon_vec = abs(epsilon_vec) - epsilon;

%% Plots theta, r and v as functions of time
figure
set(0, 'DefaultAxesFontSize', 12, 'DefaultLineLineWidth', 2)

hold on

% Plots true anomaly as a function of time
subplot(3,1,1)
plot(t_vec,theta_vec)
title('True Anomaly over Time for the Given Orbit of a Spacecraft Revolving around Earth')
xlabel('Time (s)')
ylabel('True Anomaly (radians)')

% Plots the distance from the primary body as a function of time
subplot(3,1,2)
plot(t_vec,r_vec)
title('Position over Time for the Given Orbit of a Spacecraft Revolving around Earth')
xlabel('Time (s)')
ylabel('Position (km)')

% Plots the magnitude of the velocity as a function of time
subplot(3,1,3)
plot(t_vec,v_vec)
title('Velocity over Time for the Given Orbit of a Spacecraft Revolving around Earth')
xlabel('Time (s)')
ylabel('Velocity (km/s)')

hold off

%% Calculates x- & y-distance specific values for spacecraft position at each point in time
rx = r_vec .* cos(theta_vec);
ry = r_vec .* sin(theta_vec);

% Plots the position of the spacecraft at each point in time
figure
hold on
plot(rx,ry)
title('Spacecraft Position in Orbit over Time')
xlabel('Distance (km)')
ylabel('Distance (km)')
% %% For loop defining the boundaries of Earth at each degree of true anomaly
% for theta = 1:360
%     x(theta) = 6378.1363 * cosd(theta);
%     y(theta) = 6378.1363 * sind(theta);
% end
% % Plotting Earth
% fill(x,y,'b')
hold off

%%
figure
hold on
plot(t_vec,epsilon_vec)
title('Change in Specific Energy over Time for the Given Orbit of a Spacecraft Revolving around Earth')
xlabel('Time (s)')
ylabel('Change in Specific Energy (km^2/s^2)')

hold off
