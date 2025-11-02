function [Ac, Bc, Mc] = Dyn_ibvs_c(param, hist, k)
%DYN_IBVS_C     Compute continuous-time IBVS system matrices
%
% Syntax:
%   [Ac, Bc, Mc] = Dyn_ibvs_c(param, hist, k)
%
% Description:
%   This function computes the continuous-time system matrices (Ac, Bc, Mc)
%   for an Image-Based Visual Servoing (IBVS) system at time step k.
%   It combines the image Jacobian, system dynamics, and disturbance models
%   to form the full state-space representation.
%
% Inputs:
%   param - structure containing function handles and system parameters:
%           .A6    - function handle to compute A6 matrix
%           .W     - function handle to compute disturbance matrix W
%           .Gc    - control input matrix
%           .Lsfun - function handle to compute image Jacobian Ls
%           .F1sfun- function handle to compute derivative-related matrix F1s
%   hist  - structure containing system state history:
%           .R_cl  - cell array of rotation matrices
%           .xl    - system states related to target
%           .sc    - system states related to sensor/camera
%           .w_li  - target angular velocity
%           .xs    - feature positions
%           .vc    - camera velocity/angular velocity
%           .zhat  - estimated feature depths
%   k     - current time step index
%
% Outputs:
%   Ac - continuous-time state matrix
%   Bc - continuous-time input matrix
%   Mc - continuous-time disturbance matrix
%
% Example:
%   [Ac, Bc, Mc] = Dyn_ibvs_c(param, hist, 10);
%
% Notes:
%   - Pseudo-inverse of Ls is computed with tolerance 1e-6.
%   - Image Jacobian Ls and derivative F1s are combined to form dLs.

    % Compute helper matrices for system dynamics
    A6 = param.A6(hist.R_cl{k}, hist.xl(4:6,k), hist.sc(16:18,k), hist.w_li(:,k));
    W  = param.W (hist.R_cl{k}, hist.w_li(:,k), hist.xl(7:9,k));
    Gc = param.Gc;

    % Compute image Jacobian and derivative matrices
    Ls  = param.Lsfun(hist.xs(1:8,k), hist.zhat(:,k));
    F1s = param.F1sfun(hist.xs(1:8,k), hist.vc(:,k), hist.zhat(:,k));
    dLs = F1s * Ls;

    % Compute continuous-time system matrices
    pLs = pinv(Ls, 1e-6);        % pseudo-inverse with tolerance
    Ac  = [zeros(8), eye(8); zeros(8), (dLs*pLs + Ls*A6*pLs)];
    Bc  = [zeros(8,6); Ls*Gc];
    Mc  = [zeros(8,1); Ls*W];
end
