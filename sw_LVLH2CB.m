function [s_cl, w_cl] = sw_LVLH2CB(r_t, v_t, s_ci, w_ci)
%SW_LVLH2CB     Converts the chaser's attitude and angular velocity
%               from the ECI frame to the LVLH (Local Vertical Local Horizontal)
%               frame representation.
%
% Inputs:
%   r_t  - Target position vector in the ECI frame [3×1]
%   v_t  - Target velocity vector in the ECI frame [3×1]
%   s_ci - Modified Rodrigues Parameters (MRPs) representing
%          chaser attitude from ECI to Body frame
%   w_ci - Angular velocity of the chaser expressed in the Body frame [3×1]
%
% Outputs:
%   s_cl - MRPs representing chaser attitude from LVLH to Body frame
%   w_cl - Angular velocity of the chaser relative to LVLH frame,
%          expressed in the Body frame [3×1]
%
% Description:
%   This function computes the relative attitude and angular velocity
%   of the chaser spacecraft with respect to the target’s LVLH frame.
%   The transformation is performed by combining the rotation from
%   LVLH→ECI with the chaser's attitude ECI→Body. The angular velocity
%   is adjusted by subtracting the rotation rate of the LVLH frame.


    % ---------------------------------------------------------------------
    % Step 1: Compute rotation matrices
    % ---------------------------------------------------------------------
    R_ci = mrp2dcm(s_ci);           % ECI → Body
    R_li = dcm_ECI2LVLH_rv(r_t, v_t); % ECI → LVLH
    R_cl = R_ci * R_li';            % LVLH → Body

    % ---------------------------------------------------------------------
    % Step 2: Convert relative rotation matrix to MRPs
    % ---------------------------------------------------------------------
    s_cl = dcm2mrp(R_cl);           % Body attitude relative to LVLH

    % ---------------------------------------------------------------------
    % Step 3: Compute angular velocity of the LVLH frame (in ECI)
    % ---------------------------------------------------------------------
    % The LVLH frame rotates about its z-axis with angular rate equal
    % to the orbital rate of the target satellite.
    w_li = [0; 0; norm(cross(r_t, v_t)) / norm(r_t)^2];  % [rad/s]

    % ---------------------------------------------------------------------
    % Step 4: Relative angular velocity (chaser wrt LVLH)
    % ---------------------------------------------------------------------
    w_cl = w_ci - R_cl * w_li;      % expressed in the Body frame
end
