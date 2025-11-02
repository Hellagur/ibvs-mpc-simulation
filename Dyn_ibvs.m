function param = Dyn_ibvs(param, hist, k)
%DYN_IBVS   Compute continuous and discrete IBVS system dynamic matrices
%
% Syntax:
%   param = Dyn_ibvs(param, hist, k)
%
% Description:
%   This function computes the continuous-time (Ac, Bc, Mc) and discrete-time
%   (Ad, Bd, Md) dynamic matrices for an Image-Based Visual Servoing (IBVS)
%   system at the current time step k. The discrete-time matrices are obtained
%   by discretizing the continuous-time system using the sampling period param.ts.
%
% Inputs:
%   param - structure containing system parameters, including sampling period:
%           .ts - sampling period
%   hist  - structure containing the current system state history
%   k     - current time step index
%
% Output:
%   param - updated structure with computed system matrices:
%           .Ac, .Bc, .Mc - continuous-time matrices
%           .Ad, .Bd, .Md - discrete-time matrices
%
% Example:
%   param = Dyn_ibvs(param, hist, 10);
%
% Notes:
%   - This function calls Dyn_ibvs_c to compute continuous-time matrices.
%   - Dyn_ibvs_d is used to discretize the system using param.ts.

    % Compute continuous-time system matrices
    [Ac, Bc, Mc] = Dyn_ibvs_c(param, hist, k);

    % Discretize system matrices
    [Ad, Bd, Md] = Dyn_ibvs_d(Ac, Bc, Mc, param.ts);

    % Update param structure with both continuous and discrete matrices
    param.Ac = Ac; 
    param.Bc = Bc; 
    param.Mc = Mc;
    param.Ad = Ad; 
    param.Bd = Bd; 
    param.Md = Md;
end
