


function [h_vec, r, epsilon, e_vec, e, theta, E, M_e, a, r_p, r_a, p,...
    t_p, F, M_h] = orbit_param_calc(r_vec,v_vec,mu,t)

% In case a time, t, is not input
if ~exist('t','var')
    t = 0; % time (seconds)
end

% In case mu is not input, assume mu for Earth
if ~exist('mu','var')
   mu = 398600.4415; % Standard Gravitational Parameter of Earth (km^3/s^2)
end

% Calculates the magnitue of initial vector 'r_vec'
r = sqrt(dot(r_vec,r_vec));

% Calculates specific angular momentum of initial conditions
h_vec = cross(r_vec,v_vec);

% Calculates magnitude of specific angular momentum of initial conditions
h = norm(h_vec);

% Calculates eccentricity vector for initial conditions
e_vec = (cross(v_vec,h_vec)/mu - r_vec ./ r);

% Calculates magnitude of eccentricity vector for initial conditions
e = sqrt(dot(e_vec,e_vec));

% Calculates semiparameter given initial conditions
p = dot(h_vec,h_vec) / mu;

% Calculates total orbital engergy of initial conditions
epsilon = ((dot(v_vec,v_vec) / 2) - (mu / r));

% Calculates true anomaly in radians for initial conditions
theta = acos(dot(r_vec,e_vec)/(r*e));

if e < 1
    %% Elliptical Orbits
    
    % Calculates eccentric anomaly in radians for elliptical orbits
    E = 2 * atan2(sqrt(1 - e) * tan(theta / 2), sqrt(1 + e));
    
    % Calculates mean anomaly in radians for elliptical orbits using
    % Newton-Raphson method
    M_e = E - e * sin(E);
    
    % Calculates semimajor axis length for elliptical orbits
    a = p / (1 - e^2);
    
    % Calculates length of periapsis for elliptical orbits
    r_p = a * (1 - e);
    
    % Calculates length of apoapsis for elliptical orbits
    r_a = a * (1 + e);
    
    % Calculates time since last periapsis for elliptical orbits
    t_p = t - sqrt(a^3 / mu) * M_e;
    
    % Accounting for the missing outputs if the orbit is elliptical
    M_h = 'N/A';
    F = 'N/A';
    
else
    %% Hyperbolic Orbits
    
    % Calculates eccentric anomaly in radians for hyperbolic orbit
    F = 2 * atanh(sqrt((e - 1) / (e + 1)) * tan(theta / 2));
    
    % Calculates mean anomaly in radians for hyperbolic orbit
    M_h = e * sinh(F) - F;
    
    % Calculates semimajor axis length given the initial conditions
    a = p / (e^2 - 1);
    
    % Calculates length of periapsis for initial conditions
    r_p = a * (e - 1);
    
    % Calculates length of apoapsis for initial conditions
    r_a = -a * (e + 1);
    
    % Find time since last periapsis, t_p
    t_p = h^3 / mu^2 * M_h / (e^2 - 1)^1.5;
    
    % Accounting for the missing outputs if the orbit is hyperbolic
    M_e = 'N/A';
    E = 'N/A';
end

end