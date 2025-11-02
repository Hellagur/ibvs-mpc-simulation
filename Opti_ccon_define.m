function [Aineq, bineq] = Opti_ccon_define(Ac, F, Gu, Gd, xk, delta_s, bc, mud, Rd, alpha)
%OPTI_CCON_DEFINE   Generate the deterministic equivalent of a linear chance constraint:
%                   A_c x <= b_c - b_d(k)
%                   based on the estimated disturbance statistics 
%                   (mean mu_d, covariance R_d).
%
% Inputs:
%   Ac     : Stacked state constraint matrix for the prediction horizon 
%             (size: m*Np x n*Np)
%   F      : State transition prediction matrix (n*Np x n)
%   Gu     : Control input prediction matrix (n*Np x Np*m)
%   Gd     : Disturbance propagation matrix (n*Np x Np*n)
%   xk     : Current system state (n x 1)
%   delta_s: Stacked disturbance offset vector (n*Np x 1)
%   bc     : Stacked upper bounds for state constraints (m*Np x 1)
%   mud    : Estimated mean of disturbance (n x 1)
%   Rd     : Estimated disturbance covariance (n x n)
%   alpha  : Confidence parameter (e.g., chi2inv(0.95, n))
%
% Outputs:
%   Aineq  : Equivalent deterministic inequality matrix for inputs
%   bineq  : Corresponding inequality bound adjusted by disturbance statistics
%
% The resulting constraint has the form:
%       Aineq * u_s <= bineq
% where u_s is the stacked input vector over the horizon.


    % Extract disturbance dimension and prediction horizon
    n  = size(mud, 1);             % State dimension
    m  = size(bc, 1);              % Total number of stacked inequalities
    Np = size(Gd, 2) / n;          % Prediction horizon inferred from Gd

    % Project disturbance effect on constraints
    A_Gd = Ac * Gd;                % (m x Np*n)
    b_d  = zeros(m, 1);            % Initialize tightening term

    % ---------------------------------------------------------------------
    % Compute probabilistic tightening term b_d(i)
    % Each inequality A_c(i,:) x <= b_c(i) is adjusted by the worst-case
    % standard deviation of projected disturbance:
    %
    %     b_d(i) = sum_j [ -sqrt(alpha * a_ij * R_d * a_ij') ]
    %
    % where a_ij is the part of row i of A_Gd corresponding to step j.
    % This accounts for the propagation of uncertainty at each prediction step.
    % ---------------------------------------------------------------------
    for i = 1:m
        a_i = A_Gd(i, :);
        b_sum = 0;
        for j = 1:Np
            a_ij = a_i((j-1)*n+1 : j*n);
            b_sum = b_sum - sqrt(alpha * a_ij * Rd * a_ij');  % Conservative tightening
        end
        b_d(i) = b_sum;
    end

    % ---------------------------------------------------------------------
    % Construct the deterministic linear inequality:
    %
    %   Aineq * u_s <= bineq
    %
    % where:
    %   Aineq = A_c * G_u
    %   bineq = b_c - A_c * F * x_k - A_c * G_d * δ_s + b_d
    %
    % This effectively shifts the constraint boundary inward according to
    % the estimated disturbance uncertainty.
    % ---------------------------------------------------------------------
    Aineq = Ac * Gu;
    bineq = bc - Ac * F * xk - A_Gd * delta_s + b_d;
end
