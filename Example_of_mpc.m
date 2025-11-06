clc; clear;
close all;
rng(2, "twister");  % Set random seed for reproducibility

%% ============================================================
%  Main Simulation Script for nominal IBVS-MPC
%  ------------------------------------------------------------
%  This script performs a simulation of an image-based
%  rendezvous control scenario using a nominal Model 
%  Predictive Control (MPC) framework. It initializes
%  system parameters, performs disturbance estimation,
%  solves the MPC optimization problem iteratively, and
%  visualizes results through multiple trajectory plots.
% ============================================================

%% === Step 1: Initialize system, model, and MPC parameters ===
[param, hist, mpc] = Init_params();
%  param : structure containing system constants and model parameters
%  hist  : structure storing all time histories and states
%  mpc   : structure containing MPC variables and optimization settings

%% === Step 2: Iterative MPC simulation loop ===
for k = 1:param.Tsteps
    textwaitbar(k,param.Tsteps,"MPC simulation is running");
    tcon = tic;  % Start timer for computational time measurement
    
    % ------------------------------------------------------------
    % (1) Estimate feature point depth
    % Compute the current estimated depth of image feature points
    % based on visual geometry and kinematic constraints.
    % ------------------------------------------------------------
    hist.zhat(:,k) = Est_depth(param, hist, k);

    % ------------------------------------------------------------
    % (2) Disturbance estimation
    % Update the mean and covariance of the external disturbance
    % using an exponential forgetting factor (online learning).
    % ------------------------------------------------------------
    [param.mud, param.Rd, param.gk] = Est_disturbance( ...
        hist.dT(:,k), param.mud, param.Rd, param.gk, param.lambda, k);

    % ------------------------------------------------------------
    % (3) Update system dynamics
    % Recompute Jacobians and dynamic matrices for the current step,
    % incorporating updated attitude, position, and depth estimates.
    % ------------------------------------------------------------
    param = Dyn_ibvs(param, hist, k);
    
    % ------------------------------------------------------------
    % (4) Solve the nominal MPC problem
    % Obtain the optimal control inputs and predicted trajectories
    % by minimizing the quadratic cost under deterministic constraints.
    % ------------------------------------------------------------
    [xopt, uopt] = Opti_solver(mpc, param, hist, k, @Opti_con_define_mpc);    
    hist.u(:,k) = uopt;  % Store optimal control input
    
    % ------------------------------------------------------------
    % (5) Record computation time for the MPC optimization
    % ------------------------------------------------------------
    hist.tcon(k) = toc(tcon);

    % ------------------------------------------------------------
    % (6) Propagate system dynamics
    % Simulate the next state (k+1) of the chaser spacecraft
    % using numerical integration of nonlinear dynamics.
    % ------------------------------------------------------------
    hist = Dyn_simulation(param, hist, k);
end

%% === Step 3: Post-processing and visualization ===
% The following section generates multiple trajectory plots for
% performance evaluation of the image-based MPC framework.

%% 3D relative motion in LVLH frame
Plot_relative_motion_trajectory(param, hist, k, false);

%% 2D projection (Y-X plane) of relative motion
Plot_relative_motion_trajectory2D(param, hist, k, false);

%% Feature trajectories in the camera image plane
Plot_feature_states_trajectory(param, hist, k, false);

%% Control input trajectories in the chaser body frame
Plot_control_inputs_trajectory(param, hist, k, false);

%% Relative position trajectory in the target body frame
Plot_relative_position_trajectory(param, hist, k, false);

%% Relative attitude trajectory (Modified Rodrigues Parameters)
Plot_relative_attitude_trajectory(param, hist, k, false);

%% Combined position and attitude trajectories
Plot_relative_pose_trajectories_sub(param, hist, k, false);

%% Feature error trajectories (pixel domain)
Plot_feature_error_trajectoy(param, hist, k, false);

%% Subplot view of feature tracking errors
Plot_feature_error_trajectoy_sub(param, hist, k, false);

%% Estimated depth trajectories of visual features
Plot_estimated_depth_trajectory_sub(param, hist, k, false);

%% Computation time per MPC step
Plot_time_consumption(param, hist, k, false);
