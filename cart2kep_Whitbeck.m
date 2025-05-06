function [ xKOE ] = cart2kep_Whitbeck(mu, xC)
% cart2kep : Convert inertial state vector (position and velocity) to
% Keplerian classical orbital elements.

% The input inertial state vector should be in Cartesian coordinates. For
% example, Earth-centered inertial (ECI) coordinates would be suitable to
% represent an orbit around the Earth.

%% Inertial Cartesian State Vector (6-by-1): xC
% Position
r = xC(1:3); % km
% Velocity
v = xC(4:6); % km/s

%% Vis-Viva eqn for Specific Energy
energy = dot(v,v) / 2 - mu / sqrt(dot(r,r)); % km^2/s^2

%% Specific Angular Momentum
h = cross(r,v); % km^2/s

%% Semiparameter
p = dot(h,h) / mu;

%% Eccentricity Vector
ecc = cross(v,h) / mu - r / norm(r);

%% Semi-major Axis
r_p = p / (1 + norm(ecc));
r_a = p / (1 - norm(ecc));
sma = (r_p + r_a) / 2; % km

%% Inclination
% k direction unit vector
k = [0 0 1];
% inclination
arg_inc = arccos(dot(h,k)/norm(h));
inc = acos(arg_inc); % rads

%% Line of Nodes
n = cross(k,h);
n_i = n(1);

%% True Longitude of Periapsis
arg_L = arccos(r(1) / norm(r));
L = acos(arg_L); % rads
if r(2) < 0
    L = 2*pi - L;
end

%% Argument of Latitude
arg_u = arccos(dot(n/norm(n),r) / norm(r));
u = acos(arg_u); % rads
if r(3) < 0
    u = 2*pi - u;
end

%% Longitude of Periapsis
arg_long_peri = arccos(ecc(1) / norm(ecc));
long_peri = acos(arg_long_peri);
if ecc(2) < 0
    long_peri = 2*pi - long_peri;
end

%% Handling singularities
if norm(ecc) < 10e-12 && inc < 10e-12
    RAAN = 0; % Right Ascension of Ascending Nodes (rads)
    argp = 0; % Argument of Periapsis (rads)
    tru =  L; % True Anomaly (rads)
elseif norm(ecc) < 10e-12
    argp = 0; % rads
    tru = u; % rads
elseif inc < 10e-12
    RAAN = 0; % rads
    argp = long_peri; % rads
else
    % Right Ascension of Ascending Nodes
    arg_RAAN = arccos(n_i / norm(n));
    RAAN = acos(arg_RAAN); % rads
    if n(2) < 0
        RAAN = 2*pi - RAAN; % rads
    end
    % Argument of Periapsis
    arg_argp = arccos(dot(n,ecc) / (norm(n) * norm(ecc)));
    argp = acos(arg_argp); % rads
    if ecc(3) < 0
        argp = 2*pi - argp; % rads
    end
    % True Anomaly
    tru  = acos(dot(r,ecc) / (norm(r) * norm(ecc))); % rads
    if dot(r,v) < 0
        tru = 2*pi - tru; % rads
    end
end

%% Creating xKOE Vector with Appropriate Values
xKOE = [sma; norm(ecc); inc; RAAN; argp; tru];

end
%% -----------------------------------------------------------------------+
% Author: Kendall Whitbeck
% Date: 2017-10-10
%+========================================================================+
