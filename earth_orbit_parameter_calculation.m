


function [r_mag h_vec epsilon e_vec e_mag true_anomaly_rad a r_p r_a p] = earth_orbit_parameter_calculation(r_given,v_given)

% Calculates the magnitue of initial vector 'r_vec'
r_mag = sqrt(dot(r_given,r_given))

% Calculates specific angular momentum of initial conditions
h_vec = cross(r_given,v_given)

% Calculates total orbital engergy of initial conditions
epsilon = ((dot(v_given,v_given) / 2) - (mu / r_mag))

% Calculates eccentricity vector for initial conditions
e_vec = (cross(v_given,h_vec)/mu - r_given ./ r_mag)

% Calculates magnitude of eccentricity vector for initial conditions
e_mag = sqrt(dot(e_vec,e_vec))

% Calculates true anomaly in radians for initial conditions
true_anomaly_rad = acos(dot(r_given,e_vec)/(r_mag*e_mag))

% Calculates semimajor axis length given the initial conditions
a = dot(h_vec,h_vec) / (mu * (1 - e_mag^2))

% Calculates length of periapsis for initial conditions
r_p = a * (1 - e_mag)

% Calculates length of apoapsis for initial conditions
r_a = a * (1 + e_mag)

% Calculates semiparameter given initial conditions
p = dot(h_vec,h_vec) / mu

end