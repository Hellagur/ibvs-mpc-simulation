function [traj, uopt] = Opti_solver(mpc, param, hist, k, conFun)
%OPTI_SOLVER    Solve the MPC optimization problem at time step k
%
% Inputs:
%   mpc    - structure containing MPC parameters, including prediction horizon,
%            optimization variables, objective function, and solver settings
%   param  - structure with system parameters (dynamics, state dimensions, etc.)
%   hist   - structure containing historical state and control data
%   k      - current time step index
%   conFun - function handle that builds the MPC constraint set
%
% Outputs:
%   traj - predicted optimal state trajectory (n_ftStates × (Np+1))
%   uopt - optimal control input at current step (n_controls × 1)
%
% Description:
%   This function constructs and solves a model predictive control (MPC)
%   problem using YALMIP's `optimizer` interface. It outputs the predicted
%   optimal state trajectory and control input. If the solver fails,
%   it returns the last feasible solution stored in persistent memory.
%
% Notes:
%   - The function uses persistent variables (xp, up) to retain the previous
%     feasible solution when the current optimization is infeasible.
%   - `optimizer()` allows efficient repeated MPC calls by precompiling the solver.

    persistent xp up;  % Store last feasible solution across function calls

    % Initialize outputs
    traj = zeros(param.n_ftStates, mpc.Np+1);
    uopt = zeros(param.n_controls, 1);
    
    %% --- Define MPC inputs ---
    % Typically includes current measured state and desired reference
    inputs = {hist.xs(:,k), param.sd};

    %% --- Build constraint set and cost function ---
    con = conFun(mpc, param);  % user-defined constraint builder
    
    %% --- Define MPC optimizer object ---
    % optimizer(constraints, objective, solver_options, input_vars, output_vars)
    controller = optimizer(con, mpc.obj, mpc.ops, mpc.in, mpc.out);

    %% --- Solve MPC optimization problem ---
    [sol, info] = controller(inputs);

    %% --- Check solver status ---
    if info ~= 0
        % Infeasible or failed optimization → use previous feasible result
        warning('MPC infeasible at step %d — using last feasible solution.', k);
        traj = xp;
        uopt = up;
    else
        % Successful optimization
        traj = sol{1};          % Predicted trajectory
        uopt = sol{2}(:,1);     % First control input of optimal sequence

        % Store current feasible solution for possible fallback
        xp = traj;
        up = uopt;
    end
end
