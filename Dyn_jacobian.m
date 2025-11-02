function [Lsfun, F1fun, F1sfun, F1vfun, Fifun, Fisfun] = Dyn_jacobian(param)
%DYN_JACOBIAN   Generate IBVS system Jacobian function handles
%
% Syntax:
%   [Lsfun, F1fun, F1sfun, F1vfun, Fifun, Fisfun] = Dyn_jacobian(param)
%
% Description:
%   This function generates symbolic expressions and converts them into
%   MATLAB function handles for the IBVS system Jacobians. These include:
%     - Image Jacobian (Ls)
%     - Flow-related matrix F1 and its derivatives
%     - Partial Jacobians with respect to feature positions or camera velocities
%
% Inputs:
%   param - structure containing system parameters and function handles:
%           .Lfun - function handle to compute image Jacobian for a single feature
%
% Outputs:
%   Lsfun  - function handle to compute full image Jacobian Ls(s,z)
%   F1fun  - function handle to compute flow vector F1(s,v,z)
%   F1sfun - function handle for Jacobian of F1 w.r.t feature positions s
%   F1vfun - function handle for Jacobian of F1 w.r.t camera velocity v
%   Fifun  - function handle for single-feature flow Fi(s1,s2,z1,v)
%   Fisfun - function handle for Jacobian of Fi w.r.t feature positions [s1;s2]
%
% Example:
%   [Lsfun, F1fun, F1sfun, F1vfun, Fifun, Fisfun] = Dyn_jacobian(param);
%
% Notes:
%   - Symbolic variables are used to generate the functions.
%   - Each function handle can be called with numeric inputs for efficient evaluation.
%   - This setup allows for vectorized computations of multiple features.

    % Define symbolic variables
    syms s [8 1]; syms v [6 1]; syms u [6 1]; syms z [4 1]; 
    
    % Helper for single-feature Jacobian
    L  = @(u_val, n_val, z_val) param.Lfun(u_val, n_val, z_val);

    % Compute full image Jacobian for 4 features
    Ls = [L(s(1), s(2), z(1));
          L(s(3), s(4), z(2));
          L(s(5), s(6), z(3));
          L(s(7), s(8), z(4))];

    % Flow vector
    F1 = Ls * v;

    % Convert symbolic expressions to function handles
    Lsfun  = matlabFunction(Ls, 'Vars', {s, z});
    F1fun  = matlabFunction(F1, 'Vars', {s, v, z});

    % Jacobians of F1
    F1sfun = matlabFunction(jacobian(F1, s), 'Vars', {s, v, z});
    F1vfun = matlabFunction(jacobian(F1, v), 'Vars', {s, z});

    % Single-feature flow for first feature
    Fi     = L(s(1), s(2), z(1)) * v;
    Fifun  = matlabFunction(Fi, 'Vars', {s(1), s(2), z(1), v});
    Fisfun = matlabFunction(jacobian(Fi, [s(1); s(2)]), 'Vars', {s(1), s(2), z(1), v});
end
