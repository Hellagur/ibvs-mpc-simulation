function [aDrag, tDrag] = Dyn_dis_drag(cD, m, A, we, utc, matFile, r, v, d, R_ci)
%DYN_DIS_DRAG   Compute aerodynamic drag acceleration and torque
%
% Syntax:
%   [aDrag, tDrag] = Dyn_dis_drag(cD, m, A, we, utc, matFile, r, v, d, R_ci)
%
% Description:
%   This function calculates the aerodynamic drag acceleration in the ECI frame
%   and the corresponding drag torque in the spacecraft body frame, based on
%   the satellite's position, velocity, and orientation. Atmospheric density
%   is computed using the NRLMSISE-00 model with solar and geomagnetic data.
%
% Inputs:
%   cD      - drag coefficient (dimensionless)
%   m       - satellite mass (kg)
%   A       - cross-sectional area (m^2)
%   we      - Earth's angular velocity [rad/s]
%   utc     - UTC time [year, month, day, hour, minute, second]
%   matFile - space weather data file (e.g., F10.7, geomagnetic indices)
%   r       - satellite position in ECI frame [3x1] (m)
%   v       - satellite velocity in ECI frame [3x1] (m/s)
%   d       - distance from center of mass to center of pressure (m)
%   R_ci    - rotation matrix from ECI to body frame [3x3]
%
% Outputs:
%   aDrag   - aerodynamic drag acceleration in ECI frame [3x1] (m/s^2)
%   tDrag   - aerodynamic drag torque in body frame [3x1] (N*m)
%
% Example:
%   [aDrag, tDrag] = Dyn_dis_drag(2.2, 500, 4.0, 7.2921e-5, ...
%                      [2025,11,1,12,0,0], 'spaceweather.mat', ...
%                      r, v, 0.5, R_ci);
%
% Notes:
%   - Uses eci2lla to convert ECI position to latitude, longitude, altitude.
%   - Uses fluxSolarAndGeomagnetic and atmosnrlmsise00 for atmospheric density.
%   - Drag force acts opposite to the relative velocity with respect to the rotating atmosphere.
%   - Torque is computed about the center of mass along the drag line of action.

    % --- Convert ECI to geodetic coordinates (lat, lon, alt) ---
    lla = eci2lla(r', utc);
    year = utc(1);
    dayOfYear = day(datetime(utc(1),utc(2),utc(3)), 'dayofyear');
    UTseconds = utc(4)*3600 + utc(5)*60 + utc(6);

    % --- Retrieve solar and geomagnetic indices ---
    [f107avg, f107daily, magIdx] = fluxSolarAndGeomagnetic(year, dayOfYear, UTseconds, matFile);

    % --- Compute atmospheric density using NRLMSISE-00 ---
    [~, rho] = atmosnrlmsise00(lla(3), lla(1), lla(2), year, dayOfYear, UTseconds, ...
        f107avg, f107daily, magIdx, ones(23,1), 'Oxygen', 'Warning');

    % --- Relative velocity in ECI frame (account for Earth's rotation) ---
    vrel = v - [-we*r(2); we*r(1); 0];

    % --- Aerodynamic drag force in ECI ---
    fDrag = -0.5 * cD * A * rho(6) * norm(vrel) * vrel;

    % --- Center of pressure in body frame ---
    vrel_body = R_ci * vrel;                    % ECI -> Body
    r_cp_body = -d * vrel_body / norm(vrel_body);

    % --- Transform CP to ECI frame ---
    r_cp_eci = R_ci' * r_cp_body;              % Body -> ECI

    % --- Drag torque: ECI -> Body frame ---
    tDrag_eci = cross(r_cp_eci, fDrag);
    tDrag     = R_ci * tDrag_eci;

    % --- Aerodynamic acceleration in ECI ---
    aDrag = fDrag / m;
end
