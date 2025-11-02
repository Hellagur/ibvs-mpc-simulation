function obj = Opti_obj_define(mpc)
%OPTI_OBJ_DEFINE    Construct the MPC quadratic objective function
%
% Inputs:
%   mpc - structure containing MPC configuration, including:
%         .Np : prediction horizon
%         .x  : state variables (cell array over horizon)
%         .u  : control variables (cell array over horizon)
%         .r  : reference state vector
%         .Q  : state weighting matrix
%         .R  : control weighting matrix
%         .P  : terminal weighting matrix
%
% Output:
%   obj - scalar symbolic objective function (to be minimized)
%
% Cost function form:
%   J = (x_N - r)' P (x_N - r)
%       + Σ_{k=0}^{Np-1} [ (x_k - r)' Q (x_k - r) + u_k' R u_k ]
%
% Description:
% This function defines the finite-horizon quadratic cost used in
% Model Predictive Control (MPC). The objective penalizes both
% state tracking errors and control efforts over the prediction
% horizon, with an additional terminal cost term for stability.
%
% Notes:
%   - The last term enforces a terminal penalty for asymptotic stability.
%   - An 8-dimensional zero vector is appended to r to match full state size.

    % Extract parameters
    Np = mpc.Np;
    x  = mpc.x;
    u  = mpc.u;

    % Reference trajectory (extended if necessary)
    r  = [mpc.r; zeros(8,1)];

    % Initialize objective function
    obj = 0;

    % Terminal cost
    xN = x{Np+1};
    obj = obj + (xN - r)' * mpc.P * (xN - r);

    % Stage cost (summed over prediction horizon)
    for i = 1:Np
        obj = obj + (x{i} - r)' * mpc.Q * (x{i} - r);   % state tracking term
        obj = obj + (u{i})' * mpc.R * (u{i});           % control effort term
    end
end
