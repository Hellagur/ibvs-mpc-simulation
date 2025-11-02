function con = Opti_con_define_ccmpc(mpc, param)
%OPTI_CON_DEFINE_CCMPC      Define MPC constraints with chance-constrained 
%                           state constraints in stacked form
%
% Inputs:
%   mpc   - structure containing MPC variables (x, xs, us, Np, etc.)
%   param - structure containing system matrices and constraint parameters
%
% Output:
%   con   - YALMIP-style constraint list for the optimization problem
%
% Description:
%   This function constructs the constraint set for a chance-constrained MPC
%   problem formulated in stacked (vectorized) form. 
%   It includes system dynamics, control limits, and probabilistic state constraints.

    con = [];
    Np  = mpc.Np;          % Prediction horizon
    x   = mpc.x;           % Current state variable (cell array)
    xs  = mpc.xs;          % Stacked predicted state variable
    us  = mpc.us;          % Stacked control input variable

    % Extract model parameters
    Ad  = param.Ad;        % Discrete-time system matrix A_d
    Bd  = param.Bd;        % Discrete-time input matrix B_d
    Md  = param.Md;        % Nominal disturbance vector
    mud = param.mud;       % Estimated disturbance mean
    Rd  = param.Rd;        % Estimated disturbance covariance
    alpha = param.alpha;   % Confidence level for chance constraint
    Ax  = param.Ax;        % Current-state constraint matrix
    bx  = param.bx;        % Current-state constraint bound
    Ax_ = param.Ax_;       % Stacked state constraint matrix
    bx_ = param.bx_;       % Stacked state constraint bound
    Au_ = param.Au_;       % Stacked input constraint matrix
    bu_ = param.bu_;       % Stacked input constraint bound

    % System dimensions
    nx = param.n_ftStates; % Number of state variables
    nu = param.n_controls; % Number of control inputs

    % ---------------------------------------------------------------------
    % Construct prediction matrices F, Gu, and Gd
    %   x_{k+i} = A^i x_k + Σ(A^(i-j) * B * u_{k+j-1}) + Σ(A^(i-j) * w_{k+j-1})
    % ---------------------------------------------------------------------
    F  = zeros(Np*nx, nx);     % Free response matrix
    Gu = zeros(Np*nx, Np*nu);  % Control input matrix
    Gd = zeros(Np*nx, Np*nx);  % Disturbance propagation matrix

    for i = 1:Np
        F((i-1)*nx+1:i*nx, :) = Ad^i;
        for j = 1:i
            Gu((i-1)*nx+1:i*nx, (j-1)*nu+1:j*nu) = Ad^(i-j) * Bd;
            Gd((i-1)*nx+1:i*nx, (j-1)*nx+1:j*nx) = Ad^(i-j);
        end
    end
    
    % ---------------------------------------------------------------------
    % Construct stacked disturbance sequence
    %   ds = [Md + μ_d; Md + μ_d; ...]  (Np times)
    % ---------------------------------------------------------------------
    ds = kron(ones(Np,1), Md + param.mud);

    % Get the current state (initial condition)
    xk = x{1};

    % ---------------------------------------------------------------------
    % Compute tightened chance constraints:
    % Converts probabilistic constraints P(Ax x ≤ b) ≥ α
    % into deterministic form Aineq * us ≤ bineq
    % ---------------------------------------------------------------------
    [Aineq, bineq] = Opti_ccon_define(Ax_, F, Gu, Gd, xk, ds, bx_, mud, Rd, alpha);

    % ---------------------------------------------------------------------
    % 1. Current-state constraint
    % ---------------------------------------------------------------------
    con = [con, Ax * xk <= bx];

    % ---------------------------------------------------------------------
    % 2. Dynamics constraints (stacked prediction model)
    %     xs = F*xk + Gu*us + Gd*ds
    % ---------------------------------------------------------------------
    con = [con, xs == F*xk + Gu*us + Gd*ds];

    % ---------------------------------------------------------------------
    % 3. Control input constraints
    % ---------------------------------------------------------------------
    con = [con, Au_ * us <= bu_];

    % ---------------------------------------------------------------------
    % 4. Stacked chance-constrained state inequalities
    % ---------------------------------------------------------------------
    con = [con, Aineq * us <= bineq];
end
