function con = Opti_con_define_mpc(mpc, param)
%OPTI_CON_DEFINE_MPC    Define deterministic MPC constraints in stacked 
%                       (compact) form.
%
% This function constructs the system dynamics, state, and input constraints
% over the prediction horizon, expressed in terms of stacked decision variables.
%
% The resulting constraint formulation is:
%
%   x_s = F * x_k + G_u * u_s + G_d * δ_s
%   A_u * u_s <= b_u
%   A_x * G_u * u_s <= b_x - A_x * F * x_k - A_x * G_d * δ_s
%
% Inputs:
%   mpc   : MPC structure containing symbolic optimization variables
%            - mpc.x  : current state variable
%            - mpc.xs : stacked state trajectory variable
%            - mpc.us : stacked control input variable
%            - mpc.Np : prediction horizon
%
%   param : System parameter structure
%            - Ad, Bd : discrete-time system matrices
%            - Md     : estimated disturbance mean (assumed constant)
%            - Ax, bx : single-step state constraint (at time k)
%            - Ax_, bx_: stacked state constraints (for k+1 to k+N)
%            - Au_, bu_: stacked input constraints (for k+1 to k+N)
%            - n_ftStates: number of feature states (system dimension)
%            - n_controls : number of control inputs
%
% Outputs:
%   con   : Concatenated constraint set for use in YALMIP optimizer


    % Initialize constraint set
    con = [];

    % Extract parameters
    Np  = mpc.Np;
    x   = mpc.x; 
    xs  = mpc.xs;
    us  = mpc.us;

    Ad  = param.Ad;  
    Bd  = param.Bd;  
    Md  = param.Md;

    Ax  = param.Ax;  bx  = param.bx;
    Ax_ = param.Ax_; bx_ = param.bx_;
    Au_ = param.Au_; bu_ = param.bu_;

    % System dimensions
    nx = param.n_ftStates;   % number of state variables
    nu = param.n_controls;   % number of control inputs

    % ---------------------------------------------------------------------
    % Construct stacked prediction matrices:
    %   F  : maps initial state x_k to future states
    %   Gu : maps stacked inputs to future states
    %   Gd : maps disturbances to future states
    % ---------------------------------------------------------------------
    F  = zeros(Np*nx, nx);
    Gu = zeros(Np*nx, Np*nu);
    Gd = zeros(Np*nx, Np*nx);

    for i = 1:Np
        F((i-1)*nx+1:i*nx, :) = Ad^i;
        for j = 1:i
            Gu((i-1)*nx+1:i*nx, (j-1)*nu+1:j*nu) = Ad^(i-j) * Bd;
            Gd((i-1)*nx+1:i*nx, (j-1)*nx+1:j*nx) = Ad^(i-j);
        end
    end
    
    % Stack disturbance vector: δ_s = kron(ones(Np,1), M_d)
    ds = kron(ones(Np,1), Md);

    % Get current state x_k
    xk = x{1};

    % ---------------------------------------------------------------------
    % Define MPC constraints
    % ---------------------------------------------------------------------

    % (1) Initial state constraint
    con = [con, Ax * xk <= bx];

    % (2) System dynamics constraint:
    %     x_s = F x_k + G_u u_s + G_d δ_s
    con = [con, xs == F*xk + Gu*us + Gd*ds];

    % (3) Input constraint over the horizon
    con = [con, Au_ * us <= bu_];

    % (4) Stacked state constraint:
    %     A_x * (F x_k + G_u u_s + G_d δ_s) <= b_x
    %     ⇔ A_x * G_u * u_s <= b_x - A_x * F * x_k - A_x * G_d * δ_s
    con = [con, Ax_ * Gu * us <= bx_ - Ax_ * F * xk - Ax_ * Gd * ds];
end
