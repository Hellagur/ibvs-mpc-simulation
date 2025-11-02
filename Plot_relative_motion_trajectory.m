function Plot_relative_motion_trajectory(param, hist, k, saveFig)
%PLOT_RELATIVE_MOTION_TRAJECTORY
%   Plot the relative motion trajectory of a chaser spacecraft with respect
%   to a target in LVLH (Local Vertical Local Horizontal) frame. The function
%   includes 3D visualization of:
%     - Target attitude axes
%     - Chaser trajectory and body axes
%     - Target feature points and planes
%     - Time labels
%     - Inset subplots for feature point motion and terminal state
%
% INPUTS:
%   param   - struct containing system parameters, e.g.,
%             param.ts: sampling time
%             param.xi: 3D positions of target feature points
%   hist    - struct containing simulation history, including:
%             hist.R_ti: target rotation matrices in inertial frame (cell array)
%             hist.R_li: LVLH rotation matrices (cell array)
%             hist.R_ci: chaser rotation matrices in inertial frame (cell array)
%             hist.R_cl: chaser rotation matrices in LVLH frame (cell array)
%             hist.xl  : state vector history, including chaser position
%   k       - current timestep or total number of timesteps
%   saveFig - boolean, whether to export the figure as a PDF
%
% OUTPUT:
%   None (plots are displayed and optionally saved as PDF)

    % === Define standard coordinate axes ===
    x = [1;0;0]; y = [0;1;0]; z = [0;0;1];

    % === Initialize rotation axes in LVLH frame ===
    x_lt = zeros(3, k); y_lt = zeros(3, k); z_lt = zeros(3, k);  % target axes
    x_lc = zeros(3, k); y_lc = zeros(3, k); z_lc = zeros(3, k);  % chaser axes
    for i = 1:k
        % Compute target rotation relative to LVLH
        R_tl = hist.R_ti{i} * hist.R_li{i}';
        % Compute chaser rotation relative to LVLH
        R_cl = hist.R_ci{i} * hist.R_li{i}';
        % Store rotated axes
        x_lt(:,i) = R_tl' * x;
        y_lt(:,i) = R_tl' * y;
        z_lt(:,i) = R_tl' * z;
        x_lc(:,i) = R_cl' * x;
        y_lc(:,i) = R_cl' * y;
        z_lc(:,i) = R_cl' * z;
    end

    % Chaser position history
    r_l = hist.xl(7:9,:);

    % Create main figure
    fig = figure('Units','inches','Position',[1 1 10 7]);

    %% === Main 3D Axes ===
    ax_main = axes('Position',[0.02, 0.10, 0.75, 0.85]);
    hold(ax_main, 'on'); view(ax_main, 3);
    grid(ax_main, 'on'); axis(ax_main, 'equal');

    % === Define colors ===
    color_axes  = [1 0 0; 0 1 0; 0 0 1];                % x/y/z axes
    color_point = [231 76 60; 46 204 113; 52 152 219; 241 196 15]/255; % feature points
    color_patch = [0 0 1; 1 1 0; 1 0 0];                % feature planes
    alpha_vals  = linspace(0.3, 0.5, k);                % transparency gradient

    % === Plot target axes at start and end timesteps ===
    for i = [1, k]
        start = zeros(1,3);  % LVLH origin
        Plot_arrow3D(start, start + 1.5 * x_lt(:,i)', color_axes(1,:), alpha_vals(i));
        Plot_arrow3D(start, start + 1.5 * y_lt(:,i)', color_axes(2,:), alpha_vals(i));
        Plot_arrow3D(start, start + 1.5 * z_lt(:,i)', color_axes(3,:), alpha_vals(i));
    end
    plot3(start(1),start(2),start(3),'ko','MarkerSize',3,'MarkerFaceColor','k');  % origin marker

    % === Plot chaser trajectory and body axes ===
    idx_array = [1:70:k, k];  % subsample for clarity
    for i = idx_array
        pos = r_l(:,i)';
        % Plot chaser body axes
        Plot_arrow3D(pos, pos + 1.5 * x_lc(:,i)', color_axes(1,:), alpha_vals(i));
        Plot_arrow3D(pos, pos + 1.5 * y_lc(:,i)', color_axes(2,:), alpha_vals(i));
        Plot_arrow3D(pos, pos + 1.5 * z_lc(:,i)', color_axes(3,:), alpha_vals(i));
        % Plot chaser camera model and axes
        Plot_camera_model(pos', hist.R_cl{i}', 15.0);
        Plot_camera_axes(pos', hist.R_cl{i}', 2.0, [0.2 0.6 1.0], 0.15);
    end

    % === Add time labels near cameras ===
    t_labels = [0,35,70,105,140,175,200];  
    offset_vectors = [ % offset vectors for text positioning
        1.0, -1.0,  1.0; 
       -2.0,  1.0,  1.0; 
       -1.0,  0.5,  1.5; 
       -1.0,  0.0,  1.8; 
        1.0,  0.0,  0.8; 
       -1.8,  1.0, -1.0; 
       -0.0, -1.5, -0.8  
    ];
    for j = 1:length(t_labels)
        i = t_labels(j) / param.ts + 1;
        text_pos = r_l(:,i)' + offset_vectors(j,:);
        text(ax_main, text_pos(1), text_pos(2), text_pos(3), ...
            sprintf('t = %ds', t_labels(j)), ...
            'FontSize', 10, 'FontName', 'Times New Roman', ...
            'Interpreter', 'latex', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    end

    % === Plot trajectory line ===
    plot3(r_l(1,:), r_l(2,:), r_l(3,:), 'Color', 'k', 'LineWidth', 1.5);

    %% === Plot feature planes and points ===
    xi_l = zeros(3,4);  % feature points positions
    for i = idx_array
        % Select color index
        if i == 1, idx = 1; end
        if i >= 2, idx = 2; end
        if i == k, idx = 3; end

        % Plot feature points
        R_tl = hist.R_ti{i}*hist.R_li{i}';
        for j = 1:4
            xi_l(:,j) = R_tl'*param.xi(3*(j-1)+1:3*j);
            scatter3(ax_main, xi_l(1,j), xi_l(2,j), xi_l(3,j), 30, color_point(j,:), ...
                     'filled', 'MarkerEdgeColor', color_point(j,:), ...
                     'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 1.0, 'LineWidth', 0.5);
        end

        % Plot target plane
        patch('Parent', ax_main, 'Vertices', xi_l', 'Faces', [1 2 3 4], ...
              'FaceColor', 'interp', 'FaceVertexCData', repmat(color_patch(idx,:), 4, 1), ...
              'FaceAlpha', 0.25, 'EdgeColor', 'black');
    end

    % Plot terminal light rays from chaser to feature points
    for i = 1:4
        line([r_l(1,k),xi_l(1,i)], [r_l(2,k),xi_l(2,i)], [r_l(3,k),xi_l(3,i)], ...
             'Color', color_point(i,:), 'LineStyle', '--', 'LineWidth', 0.75);
    end

    % === Axis labels ===
    xlabel('$x_L\ (\rm{m})$', 'FontSize', 12, 'FontName', 'Times New Roman', 'Interpreter', 'latex');
    ylabel('$y_L\ (\rm{m})$', 'FontSize', 12, 'FontName', 'Times New Roman', 'Interpreter', 'latex');
    zlabel('$z_L\ (\rm{m})$', 'FontSize', 12, 'FontName', 'Times New Roman', 'Interpreter', 'latex');
    set(ax_main, 'FontSize', 12, 'FontName', 'Times New Roman');

    % === Legend handles ===
    h_plane_init = patch(NaN, NaN, NaN, [0 0 1], 'FaceAlpha', 0.25, 'EdgeColor', 'black', 'DisplayName', 'Initial Feature Plane');
    h_plane_medt = patch(NaN, NaN, NaN, [1 1 0], 'FaceAlpha', 0.25, 'EdgeColor', 'black', 'DisplayName', 'Intermediate Feature Planes');
    h_plane_term = patch(NaN, NaN, NaN, [1 0 0], 'FaceAlpha', 0.25, 'EdgeColor', 'black', 'DisplayName', 'Terminal Feature Plane');
    h_traj = plot3(NaN, NaN, NaN, '-', 'Color', 'k', 'LineWidth', 1.0, 'DisplayName', 'Chaser Trajectory');
    hx = plot3(NaN, NaN, NaN, 'Color', [1 0 0], 'LineWidth', 2, 'DisplayName', '$\hat{x}$');
    hy = plot3(NaN, NaN, NaN, 'Color', [0 1 0], 'LineWidth', 2, 'DisplayName', '$\hat{y}$');
    hz = plot3(NaN, NaN, NaN, 'Color', [0 0 1], 'LineWidth', 2, 'DisplayName', '$\hat{z}$');
    legend(ax_main, [h_plane_init, h_plane_medt, h_plane_term, h_traj, hx, hy, hz], ...
           'Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 10, 'Box', 'on');

    %% === Inset subplot: Feature point motion ===
    ax_fp = axes('Position', [0.60, 0.55, 0.38, 0.38]);
    hold(ax_fp, 'on'); view(ax_fp, 3);
    plot_feature_point_motion_subplot(ax_fp, param, hist, k);

    %% === Inset: Terminal state ===
    ax_term = axes('Position', [0.60, 0.10, 0.38, 0.38]);
    plot_inset_state(ax_term, param, hist, r_l, x_lt, y_lt, z_lt, x_lc, y_lc, z_lc, k, alpha_vals, true);

    %% === Export figure as PDF if requested ===
    if saveFig == true
        set(gcf, 'PaperPositionMode', 'auto');
        set(gcf, 'Renderer', 'painters');
        figure_name = strcat('figs/relative_motion_trajectory_duration=', num2str(k*param.ts), 's');
        print(fig, figure_name, '-dpdf', '-r600')
    end
end

%% Plot target feature point motion in a subplot
function plot_feature_point_motion_subplot(ax, param, hist, k)
% plot_feature_point_motion_subplot
%   Plots the motion of the target feature points over time in a 3D inset.
%   Shows initial and terminal positions, intermediate planes, and feature point trajectories.
%
% INPUTS:
%   ax    - axes handle to draw the subplot
%   param - struct containing system parameters (param.xi: target feature points)
%   hist  - struct containing rotation matrices and state history
%   k     - total number of timesteps
%
% OUTPUT:
%   None (plots added to specified axes)

    % === Define feature point colors ===
    color_axes  = [1 0 0; 0 1 0; 0 0 1];
    color_point = [231 76 60; 46 204 113; 52 152 219; 241 196 15]/255;
    color_patch = [0 0 1; 1 1 0; 1 0 0];
    alpha_vals  = linspace(0.3, 0.5, k);

    % === Rotation axes in LVLH frame ===
    x = [1;0;0]; y = [0;1;0]; z = [0;0;1];
    x_lt = zeros(3, k); y_lt = zeros(3, k); z_lt = zeros(3, k);

    % Compute target axes in LVLH for all timesteps
    for i = 1:k
        R_tl = hist.R_ti{i} * hist.R_li{i}';
        x_lt(:,i) = R_tl' * x;
        y_lt(:,i) = R_tl' * y;
        z_lt(:,i) = R_tl' * z;
    end

    % === Plot feature planes (subset of timesteps) ===
    xi_l = zeros(3,4);
    idx_array = [1:70:k,k];
    for i = idx_array
        if i == 1, idx = 1; end
        if i >= 2, idx = 2; end
        if i == k, idx = 3; end

        % plot feature points on this plane
        R_tl = hist.R_ti{i} * hist.R_li{i}';
        for j = 1:4
            xi_l(:,j) = R_tl'*param.xi(3*(j-1)+1:3*j);
            scatter3(ax, xi_l(1,j), xi_l(2,j), xi_l(3,j), 30, color_point(j,:), 'filled', 'Marker','o', ...
                 'MarkerFaceAlpha',0.5, 'MarkerEdgeAlpha',1, 'MarkerEdgeColor', color_point(j,:), 'LineWidth',0.5);
        end
    
        % plot intermediate feature plane
        patch('Parent', ax, ...
              'Vertices', xi_l', ...
              'Faces', [1 2 3 4], ...
              'FaceColor', 'interp', ...
              'FaceVertexCData', repmat(color_patch(idx,:), 4, 1), ...
              'FaceAlpha', 0.25, 'EdgeColor', 'black');
    end

    % === Plot target axes at start and end ===
    for i = [1, k]
        start = zeros(1,3);
        Plot_arrow3D(start, start + 1.5 * x_lt(:,i)', color_axes(1,:), alpha_vals(i));
        Plot_arrow3D(start, start + 1.5 * y_lt(:,i)', color_axes(2,:), alpha_vals(i));
        Plot_arrow3D(start, start + 1.5 * z_lt(:,i)', color_axes(3,:), alpha_vals(i));
    end
    plot3(start(1), start(2), start(3), 'ko', 'MarkerSize',8, 'MarkerFaceColor','k');

    % === Set axes properties ===
    axis(ax, 'equal'); grid(ax, 'on');
    set(ax, 'FontSize', 12, 'FontName', 'Times New Roman');
end

%% Plot a small inset showing the target and chaser state
function plot_inset_state(ax, param, hist, r_l, x_lt, y_lt, z_lt, x_lc, y_lc, z_lc, idx, alpha_vals, is_terminal)
% plot_inset_state
%   Plots a 3D inset of the target and chaser state at a given timestep.
%   Shows target axes, chaser axes, feature points, and the connection lines.
%
% INPUTS:
%   ax         - axes handle for plotting the inset
%   param      - struct containing system parameters (e.g., param.xi: target feature points)
%   hist       - struct containing rotation matrices and state history
%   r_l        - chaser position history (3 x k)
%   x_lt, y_lt, z_lt - target LVLH axes (3 x k)
%   x_lc, y_lc, z_lc - chaser LVLH axes (3 x k)
%   idx        - timestep index to plot
%   alpha_vals - transparency values for target axes
%   is_terminal- boolean, true if this is a terminal state inset (for coloring)
%
% OUTPUT:
%   None (plots are added to the provided axes)

    hold(ax, 'on'); view(ax, 3); grid(ax, 'on'); axis(ax, 'equal');

    % === Define colors for axes ===
    colorX = [1 0 0]; colorY = [0 1 0]; colorZ = [0 0 1];

    % === Plot target axes at LVLH origin ===
    start = zeros(1,3);  % LVLH origin
    Plot_arrow3D(start, start + 1.5 * x_lt(:,idx)', colorX, alpha_vals(idx));
    Plot_arrow3D(start, start + 1.5 * y_lt(:,idx)', colorY, alpha_vals(idx));
    Plot_arrow3D(start, start + 1.5 * z_lt(:,idx)', colorZ, alpha_vals(idx));
    plot3(ax, start(1), start(2), start(3), 'ko', 'MarkerSize', 3, 'MarkerFaceColor', 'k');  % origin marker

    % === Plot chaser body axes at its position ===
    start = r_l(:,idx)';  % chaser position
    Plot_arrow3D(start, start + 1.5 * x_lc(:,idx)', colorX, 1.0);
    Plot_arrow3D(start, start + 1.5 * y_lc(:,idx)', colorY, 1.0);
    Plot_arrow3D(start, start + 1.5 * z_lc(:,idx)', colorZ, 1.0);
    % Plot chaser camera and its axes
    Plot_camera_model(start', hist.R_cl{idx}', 15.0);
    Plot_camera_axes(start', hist.R_cl{idx}', 2.0, [0.2 0.6 1.0], 0.15);

    % === Rotate target feature points to LVLH frame ===
    R_tl = hist.R_ti{idx} * hist.R_li{idx}';
    xi_1 = R_tl' * param.xi(1:3);
    xi_2 = R_tl' * param.xi(4:6);
    xi_3 = R_tl' * param.xi(7:9);
    xi_4 = R_tl' * param.xi(10:12);
    xi = [xi_1'; xi_2'; xi_3'; xi_4'];

    % === Plot target feature plane ===
    % Color red for terminal, blue for non-terminal
    patch('Parent', ax, 'Vertices', xi, 'Faces', [1 2 3 4], ...
          'FaceColor', 'interp', ...
          'FaceVertexCData', repmat(is_terminal * [1 0 0] + ~is_terminal * [0 0 1], 4, 1), ...
          'FaceAlpha', 0.25, 'EdgeColor', 'black');

    % === Plot feature points and connecting lines to chaser ===
    color_set = [[231 76 60]/255; [46 204 113]/255; [52 152 219]/255; [241 196 15]/255];
    marker_style = 'o';  % marker shape for feature points

    for j = 1:4
        % Plot each feature point
        scatter3(ax, xi(j,1), xi(j,2), xi(j,3), 30, color_set(j,:), 'filled', 'Marker', marker_style);
        % Draw dashed line from chaser to feature point
        line([start(1), xi(j,1)], [start(2), xi(j,2)], [start(3), xi(j,3)], ...
             'Color', color_set(j,:), 'LineStyle', '--', 'LineWidth', 0.5, 'Parent', ax);
    end

    % === Set axes properties ===
    set(ax, 'FontSize', 12, 'FontName', 'Times New Roman');
end

