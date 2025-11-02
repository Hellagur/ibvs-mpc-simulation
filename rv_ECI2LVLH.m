function [r_l, v_l] = rv_ECI2LVLH(r_t, v_t, r_c, v_c)
%RV_ECI2LVLH    Converts relative position and velocity from the ECI frame 
%               to the Local-Vertical Local-Horizontal (LVLH) frame.
%
% Inputs:
%   r_t - Target’s position vector in the ECI frame [km]
%   v_t - Target’s velocity vector in the ECI frame [km/s]
%   r_c - Chaser’s position vector in the ECI frame [km]
%   v_c - Chaser’s velocity vector in the ECI frame [km/s]
%
% Outputs:
%   r_l - Relative position vector in the LVLH frame [km]
%   v_l - Relative velocity vector in the LVLH frame [km/s]
%
% Description:
%   This function transforms the chaser’s relative motion (with respect to the target)
%   from the Earth-Centered Inertial (ECI) frame to the target’s LVLH frame.
%
%   The LVLH frame is defined as:
%       - x-axis (radial): points from Earth’s center to the target
%       - y-axis (along-track): aligned with the target’s velocity, tangent to orbit
%       - z-axis (cross-track): completes the right-hand system (normal to orbital plane)
%
%   The transformation accounts for the rotation of the LVLH frame with respect 
%   to the inertial frame using the angular velocity vector w_li = [0; 0; fDot],
%   where fDot is the target’s instantaneous angular rate along the orbit.


    % Compute DCM from ECI to LVLH frame based on target’s position and velocity
    R_li = dcm_ECI2LVLH_rv(r_t, v_t);

    % Target’s orbital angular rate (approximate instantaneous mean motion)
    fDot = norm(cross(r_t, v_t)) / norm(r_t)^2;
    w_li = [0; 0; fDot];  % Angular velocity of LVLH frame w.r.t. ECI, expressed in LVLH frame

    % Compute relative position in LVLH coordinates
    r_l  = R_li * (r_c - r_t);

    % Compute relative velocity in LVLH coordinates
    % Includes correction for frame rotation:  v_LVLH = R*v_rel_ECI - w x r
    v_l  = R_li * (v_c - v_t) - cross(w_li, r_l);
end
