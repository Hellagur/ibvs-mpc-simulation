function rhs = Dyn_ibvs_n(nmpc, param)
%DYN_IBVS_N   Compute nonlinear continuous IBVS system dynamic matrices
%
% Syntax:
%   param = Dyn_ibvs_n(param, hist, k)
%
% Description:
%   This function computes the nonlinear continuous-time dynamic matrices 
%   for an Image-Based Visual Servoing (IBVS) system.
%
% Inputs:
%   nmpc  - structure defining the NMPC optimization problem and solver settings.
%   param - structure containing system parameters.
%
% Output:
%   rhs   - nonlinear dynamics f(x,u)
%
% Example:
%   rhs = Dyn_ibvs_n(nmpc, param);


    xa = nmpc.x(1:6);       % [sigma_{CL}, omega_{CL}] in R^6
    xo = nmpc.x(7:12);      % [rho_L, rhoDot_L] in R^6
    xs = nmpc.x(13:28);     % [s, sDot] in R^16
    xz = nmpc.x(29:32);     % [depth] in R^4
    fc = nmpc.u(1:3);
    tc = nmpc.u(4:6);
    
    R_cl = mrp2dcm(xa(1:3));
    w_li = [0, 0, param.n]';
    w_ci = xa(4:6) + R_cl*w_li;

    % Compute helper matrices for feature dynamics
    A6  = param.A6(R_cl, xa(4:6), w_ci, w_li);
    W   = param.W (R_cl, w_li, xo(1:3));
    Gc  = param.Gc;
    Vc  = [R_cl*xo(4:6); xa(4:6)];

    % Compute image Jacobian and derivative matrices
    Ls  = param.Lsfun(xs(1:8), xz);
    F1s = param.F1sfun(xs(1:8), [R_cl*xo(4:6); xa(4:6)], xz);
    dLs = F1s * Ls;

    % Compute relative attitude dynamics
    rhsA = [param.Cfun(xa(1:3))*xa(4:6);
            param.A5(R_cl, w_ci, w_li)*xa(4:6) ...
            - param.Jc\param.Sfun(R_cl*w_li)*param.Jc*R_cl*w_li ...
            + param.Jc\tc];

    % Compute relative orbit dynamics
    rhsO = [xo(4:6); 
            param.A1*xo(1:3) + param.A2*xo(4:6) + R_cl'*fc];

    % Compute feature dynamics
    rhsF = [xs(9:16); 
            (Ls*A6 + dLs)*Vc + Ls*W + Ls*Gc*[fc;tc]];

    % Compute depth dynamics
    rhsD = zeros(4,1);

    rhs  = [rhsA; rhsO; rhsF; rhsD];
end
