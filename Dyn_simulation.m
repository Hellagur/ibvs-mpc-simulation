function hist = Dyn_simulation(param, hist, k)
%DYN_SIMULATION     Propagate target and chaser spacecraft states and features
%
% This function calculates the states at the next time step (k+1) for both
% target and chaser spacecraft, including relative LVLH states, feature points,
% and disturbances. It also updates the UTC time.
%
% Inputs:
%   param - structure containing spacecraft parameters, gains, and functions
%   hist  - structure storing historical states, inputs, and disturbances
%   k     - current time index
%
% Outputs:
%   hist  - updated history structure with k+1 states

    tspan = [0, param.ts];  % integration interval

    %% --- Add control input disturbances ---
    % Exponential autoregressive model for actuator errors
    Efun = @(x, alpha, sigma) alpha*x + sqrt(1-alpha^2)*sigma;
    % Rotation-based disturbance function
    Rfun = @(x) eye(3) - param.Sfun(x);
    % Apply actuator imperfections
    ufun = @(u,e,f,b,v) Rfun(e) * ((eye(3) + diag(f)) * u + b + v);
    % Saturation function
    tfun = @(u, th) min(max(u, -th*ones(3,1)), th*ones(3,1));

    % Update actuator errors for thrusters (TH) and reaction wheels (RW)
    hist.eTH(:,k+1) = Efun(hist.eTH(:,k), param.alpha_eTH, hist.eTH(:,k+1));
    hist.fTH(:,k+1) = Efun(hist.fTH(:,k), param.alpha_fTH, hist.fTH(:,k+1));
    hist.bTH(:,k+1) = Efun(hist.bTH(:,k), param.alpha_bTH, hist.bTH(:,k+1));
    hist.eRW(:,k+1) = Efun(hist.eRW(:,k), param.alpha_eRW, hist.eRW(:,k+1));
    hist.fRW(:,k+1) = Efun(hist.fRW(:,k), param.alpha_fRW, hist.fRW(:,k+1));
    hist.bRW(:,k+1) = Efun(hist.bRW(:,k), param.alpha_bRW, hist.bRW(:,k+1));

    % Compute real control inputs considering disturbances and saturations
    f_real = ufun(hist.u(1:3,k), hist.eTH(:,k+1), hist.fTH(:,k+1), hist.bTH(:,k+1), hist.vTH(:,k));
    t_real = ufun(hist.u(4:6,k), hist.eRW(:,k+1), hist.fRW(:,k+1), hist.bRW(:,k+1), hist.vRW(:,k));
    u_real = [tfun(f_real, param.fm); tfun(t_real, param.tm)];
    hist.uT(:,k) = u_real;

    %% --- Propagate spacecraft dynamics (target + chaser) ---
    % Use fixed-step RK4 integration for efficiency
    [~, x_] = odeRK4(@(t,x) Dyn_spacecrafts(t, x, u_real, param), tspan, hist.sc(:,k), 0.1);

    %% --- Extract propagated states ---
    s_ti = x_(end,1:3)';   % target attitude MRPs
    r_t  = x_(end,7:9)';   % target position ECI
    v_t  = x_(end,10:12)'; % target velocity ECI
    s_ci = x_(end,13:15)'; % chaser attitude MRPs
    w_ci = x_(end,16:18)'; % chaser angular velocity
    r_c  = x_(end,19:21)'; % chaser position ECI
    v_c  = x_(end,22:24)'; % chaser velocity ECI

    % Relative transformations: LVLH and chaser-body
    [s_cl, w_cl] = sw_LVLH2CB(r_t, v_t, s_ci, w_ci);
    [r_l, v_l]   = rv_ECI2LVLH(r_t, v_t, r_c, v_c);

    % Update UTC time
    param.utc = Dyn_updateUTC(param.utc, param.ts);

    % Store updated spacecraft states
    hist.sc(:,k+1) = x_(end,:)';
    hist.xl(:,k+1) = [s_cl; w_cl; r_l; v_l];

    %% --- Compute feature points ---
    rc = zeros(12,1);
    s  = zeros(8,1);

    % Rotation matrices
    R_ti = mrp2dcm(s_ti);
    R_ci = mrp2dcm(s_ci);
    R_cl = mrp2dcm(s_cl);
    R_ct = R_ci * R_ti';
    R_li = R_cl' * R_ci;

    hist.R_ti{k+1} = R_ti;
    hist.R_ci{k+1} = R_ci;
    hist.R_li{k+1} = R_li;
    hist.R_cl{k+1} = R_cl;
    hist.w_li(:,k+1) = [0; 0; norm(cross(r_t, v_t)) / norm(r_t)^2];

    % Compute feature depths and image-plane coordinates
    idx = @(i,n) n*(i-1)+1 : n*i;
    for j = 1:4
        rc(idx(j,3)) = param.rfun(R_cl, R_ct, r_l, param.xi(idx(j,3)));
        s(idx(j,2))  = param.Kf * rc(idx(j,3)) / rc(3*j);
    end
    vc = [R_cl*v_l; w_cl];
    ds = param.Lsfun(s, rc(3:3:end)) * vc;

    % Store feature states
    hist.xs(:,k+1) = [s; ds];
    hist.rc(:,k+1) = rc;
    hist.vc(:,k+1) = vc;

    %% --- Compute true disturbance ---
    xk1 = hist.xs(:,k+1);
    xk  = hist.xs(:,k);
    if isfield(param, 'Ad')
        dk  = xk1 - param.Ad*xk - param.Bd*hist.u(:,k) - param.Md;
        hist.dT(:,k+1) = dk;                % true disturbance
        hist.dP(:,k+1) = param.mud;         % predicted mean
        hist.dE(:,k+1) = dk - param.mud;    % estimation error
    end
end

%% ------------------------------------------------------------------------
function [t, X] = odeRK4(odefun, tspan, x0, dt)
% odeRK4  Fixed-step 4th-order Runge-Kutta integrator
%
% Usage:
%   [t, X] = odeRK4(@(t,x) f(t,x), tspan, x0, dt)
%
% Inputs:
%   odefun : function handle returning dx = f(t,x)
%   tspan  : [t0 tf] integration interval
%   x0     : initial state vector
%   dt     : integration time step
%
% Outputs:
%   t : time vector
%   X : state history (row = state at each time step)
%
% Description:
%   This simple fixed-step RK4 integrator provides fast propagation
%   for systems with moderate stiffness.

    t0 = tspan(1);
    tf = tspan(end);
    t  = (t0:dt:tf)';
    N  = length(t);
    nx = length(x0);
    X  = zeros(N, nx);
    X(1,:) = x0(:)';

    for k = 1:N-1
        tk = t(k);
        xk = X(k,:)';

        k1 = odefun(tk, xk);
        k2 = odefun(tk + 0.5*dt, xk + 0.5*dt*k1);
        k3 = odefun(tk + 0.5*dt, xk + 0.5*dt*k2);
        k4 = odefun(tk + dt, xk + dt*k3);

        X(k+1,:) = (xk + dt/6 * (k1 + 2*k2 + 2*k3 + k4))';
    end
end
