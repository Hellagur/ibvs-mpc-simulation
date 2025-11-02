function tauGrav = Dyn_dis_grav(mu, J, s, r)
%DYN_DIS_GRAV   Compute gravitational torque on a rigid satellite
%
% Syntax:
%   tauGrav = Dyn_dis_grav(mu, J, s, r)
%
% Description:
%   This function computes the gravity-gradient torque acting on a satellite
%   in orbit. The torque is calculated in the satellite body frame using
%   Modified Rodrigues Parameters (MRPs) to transform the position vector from
%   ECI to body frame.
%
% Inputs:
%   mu - gravitational constant of the central body [m^3/s^2]
%   J  - satellite inertia matrix [3x3]
%   s  - satellite's MRPs representing rotation from ECI to body frame [3x1]
%   r  - satellite position vector in ECI frame [3x1] (m)
%
% Outputs:
%   tauGrav - gravity-gradient torque in body frame [3x1] (N*m)
%
% Example:
%   tauGrav = Dyn_dis_grav(3.986e14, J, s, r);
%
% Notes:
%   - MRPs are converted to a DCM using mrp2dcm(s).
%   - Gravity-gradient torque formula: tau = 3*mu/|r|^3 * cross(r_hat, J*r_hat).

    % Transform position vector from ECI to body frame
    rhat = mrp2dcm(s) * r / norm(r);

    % Compute gravity-gradient torque
    tauGrav = 3 * mu / norm(r)^3 * cross(rhat, (J * rhat));
end
