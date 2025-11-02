function [r_c, v_c] = rv_LVLH2ECI(r_t, v_t, r_l, v_l)
%RV_LVLH2ECI    Converts relative position and velocity from the 
%               Local Vertical Local Horizontal (LVLH) frame 
%               to the Earth-Centered Inertial (ECI) frame.
%
% Syntax:
%   [r_c, v_c] = rv_LVLH2ECI(r_t, v_t, r_l, v_l)
%
% Description:
%   This function transforms the *chaser* spacecraft's position and velocity
%   from the LVLH frame (centered at the *target* spacecraft) to the ECI frame.
%   The transformation accounts for both the target’s orbital motion and the
%   rotational rate of the LVLH frame.
%
%   The relationships are given by:
%
%       r_c = r_t + R_li' * r_l
%       v_c = v_t + R_li' * ( v_l + ω_li × r_l )
%
%   where:
%       - R_li  : DCM transforming ECI → LVLH
%       - ω_li  : angular velocity of LVLH frame w.r.t. ECI, expressed in LVLH
%
% Inputs:
%   r_t : (3x1) Target position in ECI frame [m]
%   v_t : (3x1) Target velocity in ECI frame [m/s]
%   r_l : (3x1) Relative position of chaser w.r.t. target in LVLH frame [m]
%   v_l : (3x1) Relative velocity of chaser w.r.t. target in LVLH frame [m/s]
%
% Outputs:
%   r_c : (3x1) Chaser position in ECI frame [m]
%   v_c : (3x1) Chaser velocity in ECI frame [m/s]
%
% Notes:
%   - The angular rate of the LVLH frame ω_li = [0; 0; ḟ], 
%     where ḟ = |h| / |r_t|² is the orbital angular velocity magnitude.
%   - The DCM R_li is computed using the target’s orbital geometry.
%   - The returned frame transformation assumes a right-handed LVLH frame:
%         x_LVLH → radial
%         y_LVLH → along-track
%         z_LVLH → orbit normal
%
% Example:
%   r_t = [7000e3; 0; 0];
%   v_t = [0; 7.5e3; 1e3];
%   r_l = [100; 0; 0];
%   v_l = [0; 0.1; 0];
%   [r_c, v_c] = rv_LVLH2ECI(r_t, v_t, r_l, v_l)
%
% ---------------------------------------------------------------

    % --- Compute DCM from ECI → LVLH (based on target state) ---
    R_li = dcm_ECI2LVLH_rv(r_t, v_t);

    % --- Compute angular velocity of LVLH frame (about z_LVLH axis) ---
    fDot = norm(cross(r_t, v_t)) / norm(r_t)^2;  % orbital rate (rad/s)
    w_li = [0; 0; fDot];                         % ω_LVLH/ECI in LVLH frame

    % --- Transform position and velocity to ECI frame ---
    r_c = r_t + R_li' * r_l;                     % relative position
    v_c = v_t + R_li' * (v_l + cross(w_li, r_l));% relative velocity
end
