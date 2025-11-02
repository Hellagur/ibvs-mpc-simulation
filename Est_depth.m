function zhat = Est_depth(param, hist, k)
%EST_DEPTH      Estimate feature depths in IBVS using algebraic least-squares
%
% Syntax:
%   zhat = Est_depth(param, hist, k)
%
% Description:
%   This function estimates the depth (Z) of image feature points for 
%   Image-Based Visual Servoing (IBVS) using a recursive least-squares (RLS) 
%   approach on the inverse depth. If the regression vector is ill-conditioned 
%   or the estimated depth is too small, the previous estimate is retained.
%
% Inputs:
%   param - structure containing camera parameters and IBVS tuning parameters:
%           .fl  - focal length
%           .ts  - sampling period
%           .fm  - maximum allowed depth rate (used for sanity check)
%   hist  - structure containing feature history:
%           .xs    - 2n x k matrix of feature positions and velocities
%                    first n rows: positions [u1; u2; ... ; un; v1; ...]
%                    next n rows: velocities [udot1; ... ; vdotn]
%           .vc    - 6 x k camera velocity/angular velocity [vx; vy; vz; wx; wy; wz]
%           .zhat  - previous depth estimates
%   k     - current time step index
%
% Output:
%   zhat  - n x 1 vector of estimated feature depths
%
% Example:
%   zhat = Est_depth(param, hist, 10);
%
% Notes:
%   - Based on algebraic least-squares estimation of inverse depth.
%   - Incorporates threshold checks to prevent ill-conditioned estimates.


    % Extract current feature positions, velocities, and camera motion
    fl = param.fl;
    s  = hist.xs(1:8, k);    % feature positions
    ds = hist.xs(9:16, k);   % feature velocities
    vc = hist.vc(:, k);      % camera velocity/angular velocity

    n = length(s)/2;          % number of features
    zhist = 10*ones(n,1);     % default previous depth
    zhat  = zeros(n,1);       % initialize depth estimate

    if k > 1
        zhist = hist.zhat(:, k-1);
    end

    % Thresholds for regression vector and inverse depth
    tol_phi   = 1e-3;   % minimal regression vector norm squared
    tol_theta = 1e-2;   % minimal inverse-depth magnitude

    for i = 1:n
        % Parse feature inputs
        u    = s(2*(i-1)+1); udot = ds(2*(i-1)+1);
        v    = s(2*i);       vdot = ds(2*i);
        vx   = vc(1); vy = vc(2); vz = vc(3);
        wx   = vc(4); wy = vc(5); wz = vc(6);

        % Rotation-induced image flow
        Ru = (u*v/fl)*wx - ((fl^2 + u^2)/fl)*wy + v*wz;
        Rv = ((fl^2 + v^2)/fl)*wx - (u*v/fl)*wy - u*wz;

        % Residual flow
        b = [udot - Ru; vdot - Rv];

        % Regression vector for inverse depth
        phi = [-fl*vx + u*vz; -fl*vy + v*vz];

        % Check regression excitation
        if norm(phi)^2 < tol_phi
            zhat(i) = zhist(i); % retain previous estimate
            continue;
        end

        % Least-squares estimate of inverse depth
        theta = (phi' * b) / (phi' * phi);

        % Check inverse-depth magnitude
        if abs(theta) < tol_theta
            zhat(i) = zhist(i); % revert to previous
        else
            zhat(i) = 1 / theta;
        end

        % Check rate of depth change to prevent large jumps
        zrate = (zhat(i) - zhist(i)) / param.ts;
        fm = sqrt(3*param.fm^2);
        if abs(zrate) > fm*1e3
            zhat(i) = zhist(i);
        end
    end
end
