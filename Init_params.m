function [param, hist, mpc, nmpc] = Init_params(s_ct0, r_tc0, w_ti0)
%INIT_PARAMS  Initialize system, model, and MPC parameters for image-based
% rendezvous and visual servo control simulation.
%
% This function constructs the parameter structures required by the
% image-based MPC framework, including orbital dynamics, spacecraft
% properties, camera model, disturbance configuration, and optimization
% problem setup.
%
% INPUTS:
%   s_ct0  - Initial Modified Rodrigues Parameters (MRP) of the chaser
%             relative to the target body frame {T}.
%   r_tc0  - Initial relative position of the chaser in {T}-frame [m].
%   w_ti0  - (Optional) Target angular velocity in inertial frame {I} [rad/s].
%
% OUTPUTS:
%   param  - Structure containing all physical, geometric, and model parameters.
%   hist   - Structure for storing simulation histories and intermediate data.
%   mpc    - Structure defining the MPC optimization problem and solver settings.
%   nmpc   - Structure defining the NMPC optimization problem and solver settings.
%
% -------------------------------------------------------------------------


%% Target orbital parameters
% Define orbital constants and reference orbital elements
param.mu    = 398600.4418e9;                % Earth's gravitational parameter [m^3/s^2]
param.J2    = 1082.63e-6;                   % Second zonal harmonic coefficient
param.Re    = 6371e3;                       % Earth's mean radius [m]
param.we    = 7.2921e-5;                    % Earth's rotation rate [rad/s]
param.a     = param.Re + 685e3;             % Target orbit semi-major axis [m]
param.e     = 0.00155;                      % Eccentricity
param.i     = deg2rad(60);                  % Inclination [rad]
param.O     = deg2rad(30);                  % Right ascension of ascending node [rad]
param.o     = 0;                            % Argument of perigee [rad]
param.f     = 0;                            % True anomaly [rad]
param.n     = sqrt(param.mu/param.a^3);     % Mean motion [rad/s]
param.utc   = [2025, 10, 20, 12, 0, 0];     % Initial UTC time

%% Atmospheric parameters and space weather
param.cD    = 2.0;                          % Drag coefficient
param.matFile = aeroReadSpaceWeatherData('support/SW-Last5Years.csv');

%% Target spacecraft properties
param.Jt    = [17023.3,397.1,-2171.4;
                397.1,124825.7,344.2;
              -2171.4,344.2,129112.2];      % Target inertia matrix [kg·m^2]
param.mt    = 7827.8;                       % Target mass [kg]
param.Areat = 4;                            % Target cross-sectional area [m^2]
param.xi    = [0.5; 0.4;-1;
               0.5;-0.4;-1;
              -0.5;-0.4;-1;
              -0.5; 0.4;-1];                % Feature point coordinates in {T} [m]
param.d_t   = 2;                            % Distance from CoM to CP [m]

%% Chaser spacecraft properties
param.Jc    = [30.0,  0.8,  0.5;
                0.8, 32.5,  1.0;
                0.5,  1.0, 35.0];           % Chaser inertia matrix [kg·m^2]
param.mc    = 100;                          % Chaser mass [kg]
param.Areac = 1;                            % Cross-sectional area [m^2]
param.fm    = 5.0 / param.mc;               % Maximum linear acceleration [m/s^2]
param.tm    = 0.2;                          % Maximum torque [N·m]
param.ts    = 0.5;                          % Sampling period [s]
param.d_c   = 1;                            % Distance from CoM to CP [m]

% Input constraint matrices (upper/lower bounds on actuation)
param.Au    = [eye(6);-eye(6)];
param.bu    = [param.fm*ones(3,1); param.tm*ones(3,1);
               param.fm*ones(3,1); param.tm*ones(3,1)];

%% Camera and image sensor parameters
param.fl    = 20e-3/(10e-6);                % Focal length in pixel units [px]
param.f0    = 20e-3/(10e-6);                % Nominal focal length [px]
param.u0    = 0;                            % Principal point [px]
param.n0    = 0;
param.um    = 640;                          % Image half dimensions [px]
param.nm    = 512;
param.dsm   = inf;                          % No velocity constraint
param.Kf    = [param.f0, 0, param.u0; 0, param.f0, param.n0];

% State and rate constraints in image space
param.As    = [eye(8); -eye(8)];
param.bs    = [repmat([param.um; param.nm], 4, 1); repmat([param.um; param.nm], 4, 1)];
param.Ads   = [eye(8); -eye(8)];
param.bds   = param.dsm*ones(16, 1);
param.Ax    = blkdiag(param.As, param.Ads);
param.bx    = [param.bs; param.bds];

%% Helper matrix functions
% Define inline functions for skew-symmetric, attitude, and image Jacobians.
% S(v) * w = v × w skew matrix 
param.Sfun  = @(v) [  0,    -v(3),  v(2);
                    v(3),      0,  -v(1);
                   -v(2),    v(1),    0 ];
% dsdt = C(s)*w kinematics matrix of modified Rodrigues parameters
param.Cfun  = @(s) [ ...
            1 + s(1)^2 - s(2)^2 - s(3)^2,       2*(s(1)*s(2) - s(3)),           2*(s(1)*s(3) + s(2));
                2*(s(1)*s(2) + s(3)),       1 - s(1)^2 + s(2)^2 - s(3)^2,       2*(s(2)*s(3) - s(1));
                2*(s(1)*s(3) - s(2)),           2*(s(2)*s(3) + s(1)),       1 - s(1)^2 - s(2)^2 + s(3)^2 ] / 4;
% M(u,n,z) maps 3D point velocity in image space to 2×3 matrix for projection
param.Mfun  = @(u_val, n_val, z_val) [ ...
            param.f0,          0,    -(u_val - param.u0);
            0,          param.f0,    -(n_val - param.n0) ] / z_val;
% N(u,n) builds part of the image‐point dynamics (2×3)
param.Nfun  = @(u_val, n_val) [ ...
            (u_val-param.u0)*(n_val-param.n0),   - (param.f0^2 + (u_val-param.u0)^2),   param.f0*(n_val-param.n0);
            (param.f0^2 + (n_val-param.n0)^2),   - (u_val-param.u0)*(n_val-param.n0),  -param.f0*(u_val-param.u0) ] / param.f0;
% L(u,n,z) image jacobian
param.Lfun  = @(u_val, n_val, z_val) [ ...
            - param.Mfun(u_val, n_val, z_val), ...
              param.Nfun(u_val, n_val)];

% Feature-related and attitude helper functions
% r(R_cl,R_ct,r_l,xi_i) r_i = -R_cl * r_l + R_ct * xi_i
param.rfun  = @(R_cl,R_ct,r_l,xi_i) - R_cl*r_l + R_ct*xi_i;
% w(R_ct,w_tl,xi_i) w_i = R_ct * S(w_tl) * xi_i
param.wfun  = @(R_ct,w_tl,xi_i) R_ct*param.Sfun(w_tl)*xi_i;

%% Clohessy–Wiltshire (CW) dynamic matrices
% Linearized gravity gradient term
param.A1    = [3*param.n^2, 0, 0; 0, 0, 0; 0, 0, -param.n^2];
% Coriolis coupling term
param.A2    = [0, 2*param.n, 0; -2*param.n, 0, 0; 0, 0, 0];

param.A3    = @(R_cl) R_cl*param.A1*R_cl';
param.A4    = @(R_cl, w_cl) R_cl*param.A2*R_cl' - param.Sfun(w_cl);
param.A5    = @(R_cl, w_ci, w_li) -param.Jc \ param.Sfun(w_ci) * param.Jc ...
            + param.Jc \ param.Sfun(param.Jc*R_cl*w_li) - param.Sfun(R_cl*w_li);
param.A6    = @(R_cl, w_cl, w_ci, w_li) blkdiag(param.A4(R_cl, w_cl), param.A5(R_cl, w_ci, w_li));
param.A7    = @(w_cl) blkdiag(-param.S(w_cl), -param.Jc \ param.Sfun(w_cl)*param.Jc);

param.W     = @(R_cl, w_li, r_l)[param.A3(R_cl) * R_cl * r_l; 
            -param.Jc \ param.Sfun(R_cl * w_li) * param.Jc * R_cl * w_li];

param.Gc    = blkdiag(eye(3), eye(3) / param.Jc);

%% Jacobian and model functions
[Lsfun, F1fun, F1sfun, ~] = Dyn_jacobian(param);
param.Lsfun     = Lsfun; 
param.F1fun     = F1fun;
param.F1sfun    = F1sfun;

%% Initial conditions
% Compute initial orbital, attitude, and relative states of both spacecraft
% including target ECI state, LVLH transformation, and image feature states.
% Note: If target's initial angular velocity is changed, the weight matrices 
% Q and P also need to be adjusted simultaneously.

s_ti0       = [0;0;0];                      % Target's initial MRPs
if nargin < 3
    w_ti0   = deg2rad([-1.5;-1.5;1.5]);     % Target's initial angular velocity [rad/s]
end
% Target's position and velocity in {I}
[r_t0,v_t0] = rv_OE2ECI(param.mu, param.a, param.e, param.i, param.O, param.o, param.f);

s_ctd       = [0;0;0];                      % Chaser's desired relative MRPs
r_tcd       = [0;0;-5];                     % Chaser's desired relative position in {T} [m]

if nargin < 2
    s_ct0   = [0.08;0.11;0.03];             % Chaser's initial relative MRPs
    r_tc0   = [-4;3;-20];                   % Chaser's initial relative position in {T} [m]
end

R_ti0       = mrp2dcm(s_ti0);               % Initial DCM from {I} to {T}
R_li0       = dcm_ECI2LVLH_rv(r_t0, v_t0);  % Initial DCM from {I} to {L}
R_lt0       = R_li0*R_ti0';                 % Initial DCM from {T} to {L}
R_ctd       = mrp2dcm(s_ctd);               % Desired DCM from {T} to {C}
R_cld       = R_ctd*R_lt0';                 % Desired DCM from {L} to {C}
R_ct0       = mrp2dcm(s_ct0);               % Initial DCM from {T} to {C}
R_cl0       = R_ct0*R_lt0';                 % Initial DCM from {L} to {C}
R_ci0       = R_cl0 * R_li0;                % Initial DCM from {I} to {C}

r_ld        = R_lt0*r_tcd;                  % Chaser's desired relative position in {L} [m]
r_l0        = R_lt0*r_tc0;                  % Chaser's initial relative position in {L} [m]
v_l0        = [0;0;0];                      % Chaser's initial relative velocity in {L} [m/s]

s_ci0       = dcm2mrp(R_ci0);               % Chaser's initial MRPs
w_ci0       = [0; 0; 0];                    % Chaser's initial angular velocity [rad/s]
% Chaser's position and velocity in {I}
[r_c0,v_c0] = rv_LVLH2ECI(r_t0, v_t0, r_l0, v_l0);

s_cl0       = dcm2mrp(R_cl0);               % Chaser's initial relative MRPs
w_li0       = [0; 0; norm(cross(r_t0, v_t0))/norm(r_t0)^2];
w_cl0       = w_ci0 - R_cl0*w_li0;          % Chaser's initial relative angular velocity [rad/s]

rd          = zeros(12,1);                  % Feature's desired relative position in {C} [m]
rc          = zeros(12,1);                  % Feature's initial relative position in {C} [m]
idr         = @(i) 3*(i-1) + 1 : 3*i;
for i = 1:4
    rd(idr(i)) = param.rfun(R_cld,R_ctd,r_ld,param.xi(idr(i)));
    rc(idr(i)) = param.rfun(R_cl0,R_ct0,r_l0,param.xi(idr(i)));
end

sd          = zeros(8,1);                   % Initial image feature coordinates
s0          = zeros(8,1);                   % Desired image feature coordinates
ids         = @(i) 2*(i-1) + 1 : 2*i;
for i = 1:4
    sd(ids(i)) = param.Kf*rd(idr(i))/rd(3*i);
    s0(ids(i)) = param.Kf*rc(idr(i))/rc(3*i);
end
param.sd    = sd;
param.s0    = s0;

L0          = [param.Lfun(s0(1),s0(2),rc(3));
               param.Lfun(s0(3),s0(4),rc(6));
               param.Lfun(s0(5),s0(6),rc(9));
               param.Lfun(s0(7),s0(8),rc(12))];
vc          = [R_cl0*v_l0; w_cl0];          % Camera initial relative velocity
ds0         = L0*vc;                        % Initial image feature velocity

param.n_ftStates = 16;                      % The dimension of image features
param.n_scStates = 12;                      % The dimension of spacecraft
param.n_controls = 6;                       % The dimension of control inputs

%% Disturbance estimation setup
param.lambda = 0.98;                            % Forgetting factor
param.mud    = zeros(param.n_ftStates,1);       % Mean disturbance
param.Rd     = zeros(param.n_ftStates);         % Covariance
param.alpha  = chi2inv(0.95,param.n_ftStates);  % Quantiles of the chi-square distribution
param.gk     = 0;                               % Initial normalization coefficient

%% Actuator noise and bias modeling
% Define statistical parameters for reaction wheels (RW) and thrusters (TH).

% Standard deviation
param.sigma_eRW     = 1e-2;                 % Misalignment of the RWs
param.sigma_fRW     = 1e-2;                 % Scale factor bias
param.sigma_bRW     = 1e-2;                 % Offset of the RWs
param.sigma_vRW     = 1e-3;                 % Noise of the RWs
param.sigma_eTH     = 1e-3;                 % Misalignment of thrusters
param.sigma_fTH     = 1e-3;                 % Scale of the thrusters
param.sigma_bTH     = 1e-3;                 % Offset of the thrusters
param.sigma_vTH     = 1e-4;                 % Noise of the thrusters

% ECRV time constant
param.tau_eRW       = 1e8;
param.tau_fRW       = 1e8;
param.tau_bRW       = 1e8;
param.tau_eTH       = 1e8;
param.tau_fTH       = 1e8;
param.tau_bTH       = 1e8;

param.alpha_eRW     = exp(-param.ts/param.tau_eRW);
param.alpha_fRW     = exp(-param.ts/param.tau_fRW);
param.alpha_bRW     = exp(-param.ts/param.tau_bRW);
param.alpha_eTH     = exp(-param.ts/param.tau_eTH);
param.alpha_fTH     = exp(-param.ts/param.tau_fTH);
param.alpha_bTH     = exp(-param.ts/param.tau_bTH);

%% Simulation data storage
% Initialize arrays for storing state, control, and disturbance history.
param.Tsteps        = 400;

% ===== System Variables =====

% Feature's states: [s ds] in R^16 (image plane)
hist.xs             = zeros(param.n_ftStates, param.Tsteps+1);
% Target and chaser's states: 
% [sigma^T_{TI} omega^T_{TI} r^I_T v^I_T 
%  sigma^C_{CI} omega^C_{CI} r^I_C v^I_C] in R^24
hist.sc             = zeros(param.n_scStates*2, param.Tsteps+1);
% Relative states: [sigma^C_{CL} omega^C_{CL} rho^L drho^L] in R^12
hist.xl             = zeros(param.n_scStates, param.Tsteps+1);

% Control inputs: [a^C tau^C] in R^6
% T: Truth
hist.u              = zeros(param.n_controls, param.Tsteps);
hist.uT             = zeros(param.n_controls, param.Tsteps);

% ===== Auxiliary Variables =====

% Feature points position: [r1^C r2^C r3^C r4^C] in R^12
hist.rc             = zeros(12, param.Tsteps+1);
% Relative velocity: [R_{CL}*drho^L omega^C_{CL}] in R^6
hist.vc             = zeros(6, param.Tsteps+1);

% disturbace: [d_s d_ds] in R^16, k=0, 1, ...
% T: Truth, P: Predicted, E: Error
hist.dT             = zeros(param.n_ftStates, param.Tsteps+1);
hist.dP             = zeros(param.n_ftStates, param.Tsteps+1);
hist.dE             = zeros(param.n_ftStates, param.Tsteps+1);

% Rotatin matrices
hist.R_ti           = cell(1, param.Tsteps+1);
hist.R_ci           = cell(1, param.Tsteps+1);
hist.R_li           = cell(1, param.Tsteps+1);
hist.R_cl           = cell(1, param.Tsteps+1);

% Angular velocity
hist.w_li           = zeros(3, param.Tsteps+1);

% Predicted depth: [Z1^C Z2^C Z3^C Z4^C] in R^4
hist.zhat           = zeros(4, param.Tsteps+1);

% The algorithm's time consumption
hist.tcon           = zeros(1, param.Tsteps);

% Disturbance for actuators
hist.eRW            = param.sigma_eRW*randn(3, param.Tsteps+1);
hist.fRW            = param.sigma_fRW*randn(3, param.Tsteps+1);
hist.bRW            = param.sigma_bRW*randn(3, param.Tsteps+1);
hist.vRW            = param.sigma_vRW*randn(3, param.Tsteps);
hist.eTH            = param.sigma_eTH*randn(3, param.Tsteps+1);
hist.fTH            = param.sigma_fTH*randn(3, param.Tsteps+1);
hist.bTH            = param.sigma_bTH*randn(3, param.Tsteps+1);
hist.vTH            = param.sigma_vTH*randn(3, param.Tsteps);

% ===== Set Initial Values =====
hist.xs(:,1)        = [s0; ds0];
hist.sc(:,1)        = [s_ti0; w_ti0; r_t0; v_t0; s_ci0; w_ci0; r_c0; v_c0];
hist.xl(:,1)        = [s_cl0; w_cl0; r_l0; v_l0];
hist.rc(:,1)        = rc;
hist.vc(:,1)        = vc;
hist.R_ti{1}        = R_ti0;
hist.R_ci{1}        = R_ci0;
hist.R_li{1}        = R_li0;
hist.R_cl{1}        = R_cl0;
hist.w_li(:,1)      = w_li0;

%% MPC setup
% Define MPC prediction horizon, weighting matrices, and solver configuration.
mpc.Np              = 15;
mpc.Q               = diag(1e1*ones(1, param.n_ftStates));
mpc.R               = diag(1e0*ones(1, param.n_controls));
mpc.P               = diag(1e2*ones(1, param.n_ftStates));

% Define YALMIP decision variables
yalmip('clear');
nx                  = param.n_ftStates;
nu                  = param.n_controls;
Np                  = mpc.Np;

% k, ..., k+Np+1 feature states
mpc.x               = sdpvar(repmat(nx,1,Np+1),ones(1,Np+1));
% k, ..., k+Np control inputs
mpc.u               = sdpvar(repmat(nu,1,Np),ones(1,Np));
% reference position
mpc.r               = sdpvar(8,1);

% Stack vectors and define optimization objective
mpc.xs = []; for i = 2:Np+1, mpc.xs = [mpc.xs; mpc.x{i}]; end
mpc.us = []; for i = 1:Np,   mpc.us = [mpc.us; mpc.u{i}]; end

% parametersIn
mpc.in              = {mpc.x{1}, mpc.r};
% solutionsOut
mpc.out             = {[mpc.x{:}], [mpc.u{:}]};
% define objective
mpc.obj             = Opti_obj_define(mpc);

% Gurobi solver settings for numerical stability
mpc.ops = sdpsettings(...
    'solver'  , 'gurobi', ...
    'verbose' , 0, ...
    'gurobi.NumericFocus'   , 3, ...
    'gurobi.FeasibilityTol' , 1e-6, ...
    'gurobi.IntFeasTol'     , 1e-6, ...
    'gurobi.BarConvTol'     , 1e-6);

%% Constraint stacking for MPC formulation
param.Ax_   = kron(eye(Np), param.Ax);
param.bx_   = repmat(param.bx, Np, 1);
param.Au_   = kron(eye(Np), param.Au);
param.bu_   = repmat(param.bu, Np, 1);

%% NMPC setup (use CasADi/IPOPT)
import casadi.*;

nmpc.Np     = 15;
nmpc.Q      = diag(1e1*[zeros(12,1);ones(16,1);zeros(4,1)]);
nmpc.R      = diag(1e0*ones( 6,1));
nmpc.P      = diag(1e2*[zeros(12,1);ones(16,1);zeros(4,1)]);

% States: [sigma_{CL}, omega_{CL}, rho_{L}, rhoDot_{L}, s, sDot, depth]
nmpc.x      = SX.sym('x', 32, 1);
nmpc.nx     = numel(nmpc.x);

% Desired states
nmpc.Xd     = [zeros(12,1); param.sd; zeros(12,1)];

% Controls: [a^C tau^C]
nmpc.u      = SX.sym('u', 6, 1);
nmpc.nu     = numel(nmpc.u);

% Decision variables
nmpc.U      = SX.sym('U', nmpc.nu, nmpc.Np);

% Initial guess control inputs
nmpc.U0     = zeros(nmpc.Np*nmpc.nu,1);

% A matrix that represents the states over the optimization problem
nmpc.X      = SX.sym('X', nmpc.nx, nmpc.Np+1);

% Parameters (which include initial and reference states)
nmpc.X0     = SX.sym('X0', nmpc.nx*2, 1);

% Construct nonlinear dynamics
rhs         = Dyn_ibvs_n(nmpc, param);
nmpc.f      = Function('f', {nmpc.x, nmpc.u}, {rhs});

%Runge-Kutta 4
k1 = @(f,x,u,dt) ( f(x,u) );
k2 = @(f,x,u,dt) ( f(x + k1(f,x,u,dt)*dt/2,u) );
k3 = @(f,x,u,dt) ( f(x + k2(f,x,u,dt)*dt/2,u) );
k4 = @(f,x,u,dt) ( f(x + k1(f,x,u,dt)*dt,u) );
fd = @(f,x,u,dt) ( x + (dt/6) * ( k1(f,x,u,dt) + 2*k2(f,x,u,dt) + 2*k3(f,x,u,dt) + k4(f,x,u,dt) ) );

% Compute solution symbolically
nmpc.X(:,1) = nmpc.X0(1:nmpc.nx);
for k = 1:nmpc.Np
    nmpc.X(:,k+1) = fd(nmpc.f, nmpc.X(:,k), nmpc.U(:,k), param.ts);
end

% This function to get the optimal trajectory knowing the optimal solution
nmpc.F = Function('F', {nmpc.U,nmpc.X0}, {nmpc.X});

% Compute objective
nmpc.obj = 0;
r = nmpc.X0(nmpc.nx+1:nmpc.nx*2,1);
nmpc.obj = nmpc.obj + (nmpc.X(:,nmpc.Np+1) - r)'*nmpc.P*(nmpc.X(:,nmpc.Np+1) - r);
for k = 1:nmpc.Np
    nmpc.obj = nmpc.obj + (nmpc.X(:,k) - r)'*nmpc.Q*(nmpc.X(:,k) - r);
    nmpc.obj = nmpc.obj + nmpc.U(:,k)'*nmpc.R*nmpc.U(:,k);
end

% Compute constraints 
con = struct;
con.g = []; % constraint vector
for k = 1:nmpc.Np+1 % states: (s1, ..., s8)
    con.g = [con.g; nmpc.X(13:20,k)];
end

Ng = 8 * (nmpc.Np + 1);  % total number of constraints

% Preallocate
con.lbg = zeros(Ng, 1);
con.ubg = zeros(Ng, 1);

% State constraints
% For every (u_i, n_i) pair:
for k = 0:(nmpc.Np)
    idx_base = 8*k;  % base index for this time step
    % odd indices → u1,u2,u3,u4
    con.lbg(idx_base + [1 3 5 7]) = -param.um;
    con.ubg(idx_base + [1 3 5 7]) =  param.um;
    % even indices → n1,n2,n3,n4
    con.lbg(idx_base + [2 4 6 8]) = -param.nm;
    con.ubg(idx_base + [2 4 6 8]) =  param.nm;
end

% Input constraints
for k = 0:nmpc.Np-1
    idx = 6*k + (1:6);
    con.lbx(idx(1:3),1) = - param.fm;
    con.ubx(idx(1:3),1) =   param.fm;
    con.lbx(idx(4:6),1) = - param.tm;
    con.ubx(idx(4:6),1) =   param.tm;
end

nmpc.con = con;

% Set IPOPT options
opts = struct;
opts.ipopt.max_iter                     = 100;
opts.ipopt.print_level                  = 0;                % print level(0-12, 5 is verbose print)
opts.ipopt.print_timing_statistics      = 'no';             % print time statistics
opts.ipopt.print_user_options           = 'no';             % print user options
opts.print_time                         = 0;                % print solve time
opts.ipopt.output_file                  = 'log_ipopt.txt';  % output log file
opts.ipopt.acceptable_tol               = 1e-8;
opts.ipopt.acceptable_obj_change_tol    = 1e-6;

% Set IPOPT solver
n_dec       = nmpc.nu*nmpc.Np;
opt_vars    = reshape(nmpc.U, n_dec, 1);
nlp_prob    = struct('f', nmpc.obj, 'x', opt_vars, 'g', con.g, 'p', nmpc.X0);

%% ---------- sanity checks (VERY IMPORTANT) ----------
assert(length(con.lbg) == length(con.g), 'lbg length mismatch with g!');
assert(length(con.ubg) == length(con.g), 'ubg length mismatch with g!');
assert(length(con.lbx) == n_dec, 'lbx length mismatch with decision vars!');
assert(length(con.ubx) == n_dec, 'ubx length mismatch with decision vars!');

nmpc.solver = nlpsol('solver', 'ipopt', nlp_prob, opts);

end