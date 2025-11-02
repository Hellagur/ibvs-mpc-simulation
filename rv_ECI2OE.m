function oe = rv_ECI2OE(mu, r, v)
%RV_ECI2OE      Converts position and velocity vectors in the ECI frame
%               into classical orbital elements (COEs).
%
% Inputs:
%   mu - Gravitational parameter of the central body [km^3/s^2]
%   r  - Position vector of the satellite in the ECI frame [km]
%   v  - Velocity vector of the satellite in the ECI frame [km/s]
%
% Outputs:
%   oe - Vector of orbital elements [a, e, i, Ω, ω, f]:
%        a     - Semi-major axis [km]
%        e     - Eccentricity [-]
%        i     - Inclination [rad]
%        Ω     - Right ascension of ascending node (RAAN) [rad]
%        ω     - Argument of perigee [rad]
%        f     - True anomaly [rad]
%
% Description:
%   This function computes the six classical orbital elements (COEs)
%   from the satellite’s position and velocity vectors expressed in
%   the Earth-Centered Inertial (ECI) frame. The algorithm follows
%   standard orbital mechanics formulations based on conservation of
%   angular momentum and energy.
%
% Reference:
%   Vallado, D. A., *Fundamentals of Astrodynamics and Applications*, 4th ed.
%

    % ---------------------------------------------------------------------
    % Step 1: Compute specific angular momentum vector
    % ---------------------------------------------------------------------
    h = cross(r, v);                 % [km^2/s] angular momentum vector
    h_norm = norm(h);

    % ---------------------------------------------------------------------
    % Step 2: Inclination
    % ---------------------------------------------------------------------
    i = acos(h(3) / h_norm);         % inclination [rad]

    % ---------------------------------------------------------------------
    % Step 3: Right Ascension of Ascending Node (RAAN)
    % ---------------------------------------------------------------------
    K = [0; 0; 1];                   % reference z-axis of ECI frame
    n = cross(K, h);                 % node vector (points toward ascending node)
    n_norm = norm(n);
    O = atan2(n(2), n(1));           % RAAN Ω [rad]

    % ---------------------------------------------------------------------
    % Step 4: Eccentricity vector and magnitude
    % ---------------------------------------------------------------------
    e_vec = (1/mu) * ((norm(v)^2 - mu/norm(r))*r - dot(r, v)*v);
    e = norm(e_vec);                 % eccentricity magnitude

    % ---------------------------------------------------------------------
    % Step 5: Argument of perigee
    % ---------------------------------------------------------------------
    o = atan2(dot(cross(n, e_vec), h) / h_norm, dot(n, e_vec)); % ω [rad]

    % ---------------------------------------------------------------------
    % Step 6: True anomaly
    % ---------------------------------------------------------------------
    f = atan2(dot(cross(e_vec, r), h) / h_norm, dot(e_vec, r));  % f [rad]

    % ---------------------------------------------------------------------
    % Step 7: Semi-major axis
    % ---------------------------------------------------------------------
    a = h_norm^2 / mu / (1 - e^2);   % [km]

    % ---------------------------------------------------------------------
    % Output
    % ---------------------------------------------------------------------
    oe = [a; e; i; O; o; f];
end
