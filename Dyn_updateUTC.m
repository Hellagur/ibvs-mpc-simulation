function utc = Dyn_updateUTC(utc0, dt)
%DYN_UPDATEUTC      Update UTC time by a given time step
%
% Syntax:
%   utc = Dyn_updateUTC(utc0, dt)
%
% Description:
%   This function updates a UTC time vector by a specified time step (seconds),
%   using MATLAB's Julian date conversion. The output is a standard UTC vector.
%
% Inputs:
%   utc0 - initial UTC time [year, month, day, hour, minute, second]
%   dt   - time step in seconds
%
% Outputs:
%   utc  - updated UTC time [year, month, day, hour, minute, second]
%
% Example:
%   utc0 = [2025, 11, 1, 12, 0, 0];
%   dt = 10;  % 10 seconds
%   utc = Dyn_updateUTC(utc0, dt);
%
% Notes:
%   - Internally converts UTC to Julian date, adds the time step in days,
%     and converts back to date vector.

    % Convert initial UTC to Julian date
    JD = juliandate(datetime(utc0));

    % Update Julian date by dt (seconds)
    JD_new = JD + dt / 86400;  % seconds -> days

    % Convert back to UTC date vector
    utc = datevec(datetime(JD_new, 'ConvertFrom', 'juliandate'));
end
