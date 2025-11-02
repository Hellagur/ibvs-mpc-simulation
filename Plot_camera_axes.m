function Plot_camera_axes(origin, R, scale, color, alpha)
%PLOT_CAMERA_AXES   Visualize a 3D camera frustum aligned with the local z-axis.
%
% This function draws a simplified 3D camera model represented by a 
% viewing frustum and image plane, given the camera's position, orientation, 
% and scaling factor. It is useful for visualizing camera poses and 
% fields of view in spacecraft visual servoing or 3D reconstruction simulations.
%
% Inputs:
%   origin — 3×1 vector, camera center position in world coordinates.
%   R      — 3×3 rotation matrix, transformation from camera frame to world frame.
%   scale  — scalar, scaling factor defining focal length and image plane size.
%   color  — 1×3 RGB vector specifying frustum color, e.g., [0, 0, 1].
%   alpha  — scalar in [0,1], transparency of the image plane (default: 0.5).
%
% Example:
%   R = eye(3); 
%   origin = [0;0;0];
%   Plot_camera_axes(origin, R, 1.0, [0,0.7,1], 0.4);


    if nargin < 5, alpha = 0.5; end

    %% ===== Camera intrinsic geometry =====
    % Define a virtual image plane located at distance f along +z axis
    f = scale;               % focal length (depth of image plane)
    w = 0.640 * scale;       % image width
    h = 0.512 * scale;       % image height

    % Image plane corners in camera coordinates (assuming optical axis = z)
    p1 = [ w/2  h/2  f]';    % top-right corner
    p2 = [-w/2  h/2  f]';    % top-left corner
    p3 = [-w/2 -h/2  f]';    % bottom-left corner
    p4 = [ w/2 -h/2  f]';    % bottom-right corner

    % Transform image plane corners to world coordinates
    p1 = origin + R * p1;
    p2 = origin + R * p2;
    p3 = origin + R * p3;
    p4 = origin + R * p4;

    %% ===== Draw the image plane =====
    patch('XData', [p1(1), p2(1), p3(1), p4(1)], ...
          'YData', [p1(2), p2(2), p3(2), p4(2)], ...
          'ZData', [p1(3), p2(3), p3(3), p4(3)], ...
          'FaceColor', color, 'FaceAlpha', alpha, ...
          'EdgeColor', 'k');  % black edge for better visibility

    %% ===== Draw frustum lines =====
    % Connect camera center with image plane corners
    line([origin(1) p1(1)], [origin(2) p1(2)], [origin(3) p1(3)], 'Color', color);
    line([origin(1) p2(1)], [origin(2) p2(2)], [origin(3) p2(3)], 'Color', color);
    line([origin(1) p3(1)], [origin(2) p3(2)], [origin(3) p3(3)], 'Color', color);
    line([origin(1) p4(1)], [origin(2) p4(2)], [origin(3) p4(3)], 'Color', color);

    %% ===== Draw camera center =====
    % Represent the camera optical center as a small black sphere
    plot3(origin(1), origin(2), origin(3), 'ko', ...
          'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'MarkerSize', 2);
end
