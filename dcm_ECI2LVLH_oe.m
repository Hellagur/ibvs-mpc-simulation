function R_li = dcm_ECI2LVLH_oe(i, O, o, f)
%DCM_ECI2LVLH_OE    Direction cosine matrix from ECI to LVLH frame
%
% Syntax:
%   R_li = dcm_ECI2LVLH_oe(i, O, o, f)
%
% Description:
%   Computes the direction cosine matrix (DCM) transforming a vector from
%   the Earth-Centered Inertial (ECI) frame to the Local Vertical Local
%   Horizontal (LVLH) frame, using classical orbital elements:
%       - inclination (i)
%       - right ascension of ascending node (O)
%       - argument of perigee (o)
%       - true anomaly (f)
%
% Inputs:
%   i  - inclination [rad]
%   O  - right ascension of ascending node (RAAN) [rad]
%   o  - argument of perigee [rad]
%   f  - true anomaly [rad]
%
% Outputs:
%   R_li - DCM from ECI to LVLH frame (3×3)
%
% Notes:
%   The LVLH frame is defined as:
%       x_LVLH : toward the velocity direction
%       z_LVLH : toward the radial direction (nadir)
%       y_LVLH : completes right-hand triad
%
%   The ECI to rotation from orbital plane is given by:
%       R_LVLH_ECI = R3(O) * R1(i) * R3(o + f)
%
% Example:
%   i = deg2rad(51.6); O = deg2rad(120); o = deg2rad(40); f = deg2rad(10);
%   R_li = dcm_ECI2LVLH_oe(i, O, o, f);
%
% See also: angle2dcm

    % Rotation from orbital frame to ECI
    R_li = angle2dcm(O, i, o + f, 'ZXZ');
end
