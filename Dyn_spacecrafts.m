function dx = Dyn_spacecrafts(~, x, u, param)
%DYN_SPACECRAFTS    Compute satellite attitude and orbit dynamics
%
% Syntax:
%   dx = Dyn_spacecrafts(t, x, u, param)
%
% Description:
%   This function computes the derivatives of states for a target and chaser
%   satellite pair, including attitude (MRPs), angular velocity, position,
%   and velocity. It accounts for gravity-gradient torque, J2 perturbation,
%   and aerodynamic drag.
%
% Inputs:
%   ~     - time (not used, included for ODE solver compatibility)
%   x     - state vector [24x1] containing target and chaser states:
%           target: st(1:3), wt(4:6), rt(7:9), vt(10:12)
%           chaser: sc(13:15), wc(16:18), rc(19:21), vc(22:24)
%   u     - control input vector [6x1], e.g., u(1:3) = translational, u(4:6) = rotational
%   param - structure containing physical constants and function handles:
%           .mu      - gravitational constant
%           .J2      - Earth's J2 coefficient
%           .Re      - Earth's equatorial radius
%           .we      - Earth's angular velocity
%           .utc     - UTC time vector [year,month,day,hour,minute,second]
%           .matFile - space weather data file
%           .Cfun    - function handle: MRPs to DCM derivative matrix
%           .Sfun    - function handle: skew-symmetric matrix of angular velocity
%           .Jt, .mt, .Areat - target inertia, mass, area
%           .Jc, .mc, .Areac - chaser inertia, mass, area
%           .d_t, .d_c       - center-of-pressure offset for target/chaser
%
% Outputs:
%   dx - derivative of state vector [24x1]:
%        [st_dot; wt_dot; rt_dot; vt_dot; sc_dot; wc_dot; rc_dot; vc_dot]
%
% Example:
%   dx = Dyn_spacecrafts(0, x0, u, param);
%
% Notes:
%   - Uses MRPs for attitude representation.
%   - Gravity-gradient torque, J2 acceleration, and aerodynamic drag are computed
%     for both target and chaser satellites.
%   - Control inputs are applied to the chaser only.

    %% Extract constants and parameters
    mu      = param.mu;
    J2      = param.J2;
    cD      = param.cD;
    Re      = param.Re;
    we      = param.we;
    utc     = param.utc;
    matFile = param.matFile;
    Cfun    = param.Cfun;
    Sfun    = param.Sfun;

    Jt      = param.Jt;
    mt      = param.mt;
    At      = param.Areat;

    Jc      = param.Jc;
    mc      = param.mc;
    Ac      = param.Areac;

    %% Extract target states
    st = x(1:3);  wt = x(4:6);
    rt = x(7:9);  vt = x(10:12);

    %% Compute target dynamics
    R_ti    = mrp2dcm(st);
    tGrav_t = Dyn_dis_grav(mu, Jt, st, rt);
    aJ2_t   = Dyn_dis_J2(mu, J2, Re, rt);
    [aDrag_t, tDrag_t] = Dyn_dis_drag(cD, mt, At, we, utc, matFile, rt, vt, param.d_t, R_ti);

    %% Extract chaser states
    sc = x(13:15); wc = x(16:18);
    rc = x(19:21); vc = x(22:24);

    %% Compute chaser dynamics
    R_ci    = mrp2dcm(sc);
    tGrav_c = Dyn_dis_grav(mu, Jc, sc, rc);
    aJ2_c   = Dyn_dis_J2(mu, J2, Re, rc);
    [aDrag_c, tDrag_c] = Dyn_dis_drag(cD, mc, Ac, we, utc, matFile, rc, vc, param.d_c, R_ci);

    %% State derivatives
    dx = zeros(24,1);

    % Target attitude and angular velocity
    dx(1:3)   = Cfun(st) * wt;
    dx(4:6)   = Jt \ (-Sfun(wt)*Jt*wt + tGrav_t + tDrag_t);

    % Target position and velocity
    dx(7:9)   = vt;
    dx(10:12) = -mu*rt/norm(rt)^3 + aDrag_t + aJ2_t;

    % Chaser attitude and angular velocity
    dx(13:15) = Cfun(sc) * wc;
    dx(16:18) = Jc \ (-Sfun(wc)*Jc*wc + u(4:6) + tGrav_c + tDrag_c);

    % Chaser position and velocity
    dx(19:21) = vc;
    dx(22:24) = -mu*rc/norm(rc)^3 + R_ci' * u(1:3) + aDrag_c + aJ2_c;
end
