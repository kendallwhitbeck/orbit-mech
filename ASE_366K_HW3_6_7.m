% clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 6 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Given:
mu = 1;
r_given = [-1.1545, 1.2795, -.3859];
v_given = [-.52992, -.22663, .4896];

% Calculates the magnitue of initial vector 'r_vec'
r_mag = sqrt(dot(r_given,r_given));

% Calculates specific angular momentum of initial conditions
h_vec = cross(r_given,v_given);

% Calculates total orbital engergy of initial conditions
epsilon = ((dot(v_given,v_given) / 2) - (mu / r_mag));

% Calculates eccentricity vector for initial conditions
e_vec = (cross(v_given,h_vec)/mu - r_given ./ r_mag);

% Calculates magnitude of eccentricity vector for initial conditions
e_mag = sqrt(dot(e_vec,e_vec));

% Calculates true anomaly in radians for initial conditions
theta_rad = acos(dot(r_given,e_vec)/(r_mag*e_mag));

% Calculates semi-major axis length given the initial conditions
a = dot(h_vec,h_vec) / (mu * (1 - e_mag^2));

% Calculates length of periapsis for initial conditions
r_p = a * (1 - e_mag);

% Calculates length of apoapsis for initial conditions
r_a = a * (1 + e_mag);

% Calculates semiparameter given initial conditions
p = dot(h_vec,h_vec) / mu;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 7 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% True anomaly vector every degree from 1 to 360
true_anomaly = 1:360;

% Orbit equation as a function of true anomaly
r_vec = p ./ (1 + e_mag .* cos(true_anomaly.*pi./180));

% Calculates orbital velocity for each degree of true anomaly
v_vec = sqrt(mu * (2 ./ r_vec) - (1 / a));

% Calculates specific kinetic energy for each degree of true anomaly
spec_kinetic_E_vec = v_vec .^ 2 ./ 2;

% Calculates specific potential energy for each degree of true anomaly
spec_potential_E_vec = mu ./ r_vec;

% Calculates total specific energy for each degree of true anomaly
total_spec_E_vec = spec_kinetic_E_vec - spec_potential_E_vec;

% Plots specific kinetic energy, specific potential energy and total
% specific energy
figure
hold on
plot(true_anomaly,spec_kinetic_E_vec)
plot(true_anomaly,spec_potential_E_vec)
plot(true_anomaly,total_spec_E_vec, 'g')
set(0, 'DefaultAxesFontSize', 12, 'DefaultLineLineWidth', 3)
title('The Kinetic, Potential & Total Energies for an Elliptical Orbit defined by Semiparameter p =  1.7660 (DU)')
xlabel('True Anomaly (degrees)')
ylabel('Specific Energy (DU^2/TU^2)')
xlim([0 360])
ylim([-.35 .65])
legend('Specific Kinetic Energy','Specific Potential Energy','Total Specific Energy')
hold off

% My results agree with what I expect: the kinetic and potential specific
% energies vary sinusoidally depending on the satellite's proximity to the
% primary body and the total specific energy is negative & remains constant
% throughout the entire orbit. Kepler's 2nd Law states that equal area must
% be covered over equal time, implying that satellites move faster
% when close to the primary body. This agrees with the fact that the
% kinetic energy increase as the satellite approaches periapsis. The
% increase in potential energy as periapsis is neared is due to the inverse
% relationship of potential energy and distance from the primary body. The
% total specific energy remains constant as predicted by the law of the 
% conservation of energy. Lastly, the total specific energy is negative, as
% should be the case for a satellite trapped in an elliptical or circular
% orbit.




