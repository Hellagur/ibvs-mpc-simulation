function [mu_d, R_d, g_k] = Est_disturbance(Dk_hat, mu_prev, R_prev, g_prev, lambda, k)
%EST_DISTURBANCE    Recursively estimates Gaussian disturbance mean and covariance
%
% Syntax:
%   [mu_d, R_d, g_k] = Est_disturbance(Dk_hat, mu_prev, R_prev, g_prev, lambda, k)
%
% Description:
%   This function implements a recursive estimation of the disturbance mean
%   and covariance based. A forgetting factor lambda is applied to weigh 
%   recent data more heavily.
%
% Inputs:
%   Dk_hat  - (n x 1 double) Current disturbance residual
%   mu_prev - (n x 1 double) Previous disturbance mean estimate
%   R_prev  - (n x n double) Previous disturbance covariance estimate
%   g_prev  - (1 x 1 double) Previous normalization factor
%   lambda  - (1 x 1 double) Forgetting factor (0 < lambda <= 1)
%   k       - (1 x 1 double) Current time step index
%
% Outputs:
%   mu_d    - (n x 1 double) Current disturbance mean estimate
%   R_d     - (n x n double) Current disturbance covariance estimate
%   g_k     - (1 x 1 double) Updated normalization factor
%
% Example:
%   [mu_d, R_d, g_k] = Est_disturbance(Dk_hat, mu_prev, R_prev, g_prev, 0.98, 10);


% Update normalization factor
g_k = g_prev + exp(-lambda * k);

% Compute exponential weighting term
expT = exp(-lambda) / g_k;

% Update mean estimate
mu_d = expT * (g_prev * mu_prev + Dk_hat);

% Compute deviation from updated mean
delta = Dk_hat - mu_d;

% Update covariance estimate
R_d = expT * (g_prev * R_prev + delta * delta');
end