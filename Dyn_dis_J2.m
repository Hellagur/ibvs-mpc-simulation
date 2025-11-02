function aJ2 = Dyn_dis_J2(mu, J2, Re, r)
%DYN_DIS_J2     Compute J2 perturbation acceleration in ECI frame
%
% Syntax:
%   aJ2 = Dyn_dis_J2(mu, J2, Re, r)
%
% Description:
%   This function computes the acceleration due to the Earth's second zonal
%   harmonic (J2) perturbation for a satellite at position r in the ECI frame.
%   The effect of the equatorial bulge is included in the acceleration vector.
%
% Inputs:
%   mu  - gravitational constant of Earth [m^3/s^2]
%   J2  - Earth's second zonal harmonic coefficient (dimensionless)
%   Re  - Earth's equatorial radius [m]
%   r   - satellite position in ECI frame [3x1] (m)
%
% Outputs:
%   aJ2 - J2 perturbation acceleration in ECI frame [3x1] (m/s^2)
%
% Example:
%   aJ2 = Dyn_dis_J2(3.986e14, 1.08263e-3, 6378137, r);
%
% Notes:
%   - The acceleration formula is based on standard gravity model including J2:
%       aJ2 = -1.5 * J2 * mu * Re^2 / |r|^5 * [x*(1-5*z^2/|r|^2); 
%                                              y*(1-5*z^2/|r|^2);
%                                              z*(3-5*z^2/|r|^2)]
%   - r = [x; y; z] in ECI frame.

    % Precompute powers and ratios
    rNorm = norm(r);
    z_r2 = (r(3)/rNorm)^2;

    % Compute J2 acceleration
    aJ2 = -1.5 * J2 * mu * Re^2 / rNorm^5 * ...
           [r(1)*(1 - 5*z_r2);
            r(2)*(1 - 5*z_r2);
            r(3)*(3 - 5*z_r2)];
end
