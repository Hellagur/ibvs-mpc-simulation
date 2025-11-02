function [r, v] = rv_OE2ECI(mu, a, e, i, O, o, f)
%RV_OE2ECI      Converts classical orbital elements (COEs) to position and velocity 
%               vectors in the Earth-Centered Inertial (ECI) frame.
%
% Inputs:
%   mu - Gravitational parameter [km^3/s^2]
%   a  - Semi-major axis [km]
%   e  - Eccentricity [-]
%   i  - Inclination [rad]
%   O  - Right ascension of ascending node (RAAN) [rad]
%   o  - Argument of perigee [rad]
%   f  - True anomaly [rad]
%
% Outputs:
%   r  - Position vector in the ECI frame [km]
%   v  - Velocity vector in the ECI frame [km/s]


    % ----------------------------
    % Compute orbital radius:
    % r_p = a(1 - e^2) / (1 + e*cos(f))
    % ----------------------------
    r = a*(1 - e^2) / (1 + e*cos(f)) * ...
        [cos(f+o)*cos(O) - cos(i)*sin(f+o)*sin(O);
         cos(i)*cos(O)*sin(f+o) + cos(f+o)*sin(O);
         sin(i)*sin(f + o)];
    % Result: r is the position vector in ECI coordinates

    % ----------------------------
    % Compute orbital velocity vector:
    % v = sqrt(mu/p) * [ ... ]
    % where p = a(1 - e^2)
    % This formulation directly transforms velocity from the orbital plane to ECI.
    % ----------------------------
    v = sqrt(mu / (a*(1 - e^2))) * ...
        [-cos(O)*sin(f+o) - sin(O)*cos(i)*cos(f+o) - e*(cos(O)*sin(o) + sin(O)*cos(o)*cos(i));
         -sin(O)*sin(f+o) + cos(O)*cos(i)*cos(f+o) - e*(sin(O)*sin(o) - cos(O)*cos(o)*cos(i));
          sin(i) * (cos(f+o) + e*cos(o))];
    % Result: v is the velocity vector in ECI coordinates
end
