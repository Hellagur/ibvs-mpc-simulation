function [Ad, Bd, Md] = Dyn_ibvs_d(Ac, Bc, Mc, ts)
%DYN_IBVS_D     Discretize continuous-time IBVS system matrices
%
% Syntax:
%   [Ad, Bd, Md] = Dyn_ibvs_d(Ac, Bc, Mc, ts)
%
% Description:
%   This function computes the discrete-time system matrices (Ad, Bd, Md)
%   for an Image-Based Visual Servoing (IBVS) system by applying the
%   state-transition matrix method (matrix exponential) to the augmented
%   continuous-time system [Ac, Bc, Mc].
%
% Inputs:
%   Ac - continuous-time state matrix (n x n)
%   Bc - continuous-time input matrix (n x m)
%   Mc - continuous-time disturbance matrix (n x 1)
%   ts - sampling period
%
% Outputs:
%   Ad - discrete-time state matrix (n x n)
%   Bd - discrete-time input matrix (n x m)
%   Md - discrete-time disturbance matrix (n x 1)
%
% Example:
%   [Ad, Bd, Md] = Dyn_ibvs_d(Ac, Bc, Mc, 0.01);
%
% Notes:
%   - The augmented matrix [Ac, Bc, Mc; zeros(m+1, n+m+1)] is exponentiated
%     to compute the discrete-time transition.
%   - Md corresponds to the last column of the augmented matrix exponential.

    % Dimensions of the continuous-time system
    n = size(Ac, 1);
    m = size(Bc, 2);

    % Construct augmented matrix for state + input + disturbance
    M = [Ac, Bc, Mc;
         zeros(m+1, n+m+1)];

    % Compute matrix exponential for discretization
    phi = expm(M * ts);

    % Extract discrete-time matrices
    Ad = phi(1:n, 1:n);
    Bd = phi(1:n, n+1:n+m);
    Md = phi(1:n, end);
end
