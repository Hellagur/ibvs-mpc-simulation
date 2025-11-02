function R_li = dcm_ECI2LVLH_rv(r, v)
%DCM_ECI2LVLH_RV    Computes the direction cosine matrix (DCM) 
%                   from the Earth-Centered Inertial (ECI) frame 
%                   to the Local Vertical Local Horizontal (LVLH) frame.
%
% Syntax:
%   R_li = dcm_ECI2LVLH_rv(r, v)
%
% Description:
%   The LVLH frame (also known as the Hill frame or orbit frame) is defined as:
%       - x_LVLH : Radial direction, pointing from the Earth to the spacecraft.
%       - y_LVLH : Along-track direction, tangential to the orbit and aligned with
%                  the spacecraft's instantaneous velocity vector (prograde).
%       - z_LVLH : Orbit normal direction, completing the right-handed triad.
%
%   The returned rotation matrix R_li transforms a vector from the ECI frame 
%   to the LVLH frame, such that:
%       v_LVLH = R_li * v_ECI
%
% Inputs:
%   r : (3x1) Position vector of the spacecraft in the ECI frame [m]
%   v : (3x1) Velocity vector of the spacecraft in the ECI frame [m/s]
%
% Outputs:
%   R_li : (3x3) Direction Cosine Matrix transforming ECI -> LVLH
%
% Notes:
%   - The resulting frame is right-handed, satisfying:
%         i_r × i_theta = i_h
%   - Numerical normalization is included for robustness against 
%     floating-point errors in cross products.
%   - Some definitions use the nadir direction (-r̂) for x_LVLH; modify the 
%     sign of 'ir' accordingly if a nadir-pointing frame is desired.
%
% Example:
%   r = [7000e3; 0; 0]; 
%   v = [0; 7.5e3; 1.0];
%   R_li = dcm_ECI2LVLH_rv(r, v)
%
%   % Verify orthogonality:
%   disp(R_li * R_li');  % Should be approximately identity
% ---------------------------------------------------------------

    % --- Compute orbital angular momentum (normal to the orbit plane) ---
    h  = cross(r, v);          % [m^2/s], points along orbit normal

    % --- Define unit vectors of LVLH frame in ECI coordinates ---
    ir = r / norm(r);          % x-axis: radial direction (local vertical)
    ih = h / norm(h);          % z-axis: orbit normal (local normal)
    itheta = cross(ih, ir);    % y-axis: along-track direction (local horizontal)
    itheta = itheta / norm(itheta);  % ensure orthogonality numerically

    % --- Assemble DCM (rows are LVLH unit vectors in ECI coordinates) ---
    R_li = [ir, itheta, ih]';  % transpose to get ECI -> LVLH mapping
end
