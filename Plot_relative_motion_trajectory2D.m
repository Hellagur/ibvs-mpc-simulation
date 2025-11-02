function Plot_relative_motion_trajectory2D(param, hist, k, saveFig)
% Plot_relative_motion_trajectory
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
    fig = figure('Units','inches','Position',[1 1 8 6]);

    %% === Main 3D Axes ===
    ax_main = axes('Position',[0.05, 0.05, 0.90, 0.90]);
    hold(ax_main, 'on'); view(ax_main, -90, 90);
    grid(ax_main, 'on'); axis(ax_main, 'equal');
    xlim([-3.5,7.5]); ylim([-16.5,5.5])

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
    idx_array = [1,71,141,k];  % subsample for clarity
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
    t_labels = [0, 35, 70, 200];  
    offset_vectors = [ % offset vectors for text positioning
        1.5,  0.2, 0.0; 
       -1.2, -0.5, 0.0; 
        1.5, -0.0, 0.0; 
       -0.2, -1.6, 0.0; 
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

        % Plot light rays from chaser to feature points
        for j = 1:4
            line([r_l(1,i),xi_l(1,j)], [r_l(2,i),xi_l(2,j)], [r_l(3,i),xi_l(3,j)], ...
                 'Color', color_point(j,:), 'LineStyle', '--', 'LineWidth', 0.75);
        end
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

    %% === Export figure as PDF if requested ===
    if saveFig == true
        set(gcf, 'PaperPositionMode', 'auto');
        set(gcf, 'Renderer', 'painters');
        figure_name = strcat('figs/relative_motion_trajectory2D_duration=', num2str(k*param.ts), 's');
        print(fig, figure_name, '-dpdf', '-r600')
    end
end