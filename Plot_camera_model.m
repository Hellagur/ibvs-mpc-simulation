function Plot_camera_model(origin, R, scale, bodyColor, lensColor)
%PLOT_CAMERA_MODEL  Draw a 3D camera model (body + lens) at a specified pose.
%
% This function visualizes a 3D camera body composed of a rectangular 
% cube (main body) and a short cylinder (lens). The model is oriented 
% along the camera’s local +z axis and transformed to world coordinates 
% according to the given rotation matrix and position vector.
%
% Inputs:
%   origin     — 3×1 vector, camera center position in world coordinates.
%   R          — 3×3 rotation matrix, camera-to-world transformation.
%   scale      — scalar, scaling factor controlling overall model size (default: 1.0).
%   bodyColor  — 1×3 RGB vector specifying camera body color (default: [0.8 0.8 0.8]).
%   lensColor  — 1×3 RGB vector specifying lens color (default: [0.5 0.5 0.5]).
%
% Example:
%   origin = [0;0;0];
%   R = eye(3);
%   Plot_camera_model(origin, R, 1.0, [0.7 0.7 0.7], [0.2 0.2 0.8]);
%
% Description:
%   - The main body is represented by a rectangular cube centered at
%     (0, 0, -0.025) in the camera frame.
%   - The lens is modeled as a short cylinder aligned with the +z axis.
%   - The combined model provides an intuitive visualization of 
%     camera pose and orientation in 3D scenes, e.g., in spacecraft
%     visual servoing or SLAM simulation.

    %% ===== Default Parameters =====
    if nargin < 5
        lensColor = [0.5 0.5 0.5];
    end
    if nargin < 4
        bodyColor = [0.8 0.8 0.8];
    end
    if nargin < 3
        scale = 1.0;
    end

    %% ===== Create Camera Body (Cube) =====
    % Cube centered around (0,0,-0.025) in the camera frame
    [Xb, Yb, Zb] = meshgrid([-0.025 0.025], [-0.025 0.025], [-0.03 0]);
    Xb = scale * Xb(:)'; 
    Yb = scale * Yb(:)';
    Zb = scale * Zb(:)';

    % Transform vertices into world coordinates
    verts = [Xb; Yb; Zb];
    verts = R * verts + origin;

    % Define cube faces
    faces = [
        1 2 4 3;   % Front
        5 6 8 7;   % Back
        1 2 6 5;   % Top
        3 4 8 7;   % Bottom
        1 3 7 5;   % Left
        2 4 8 6    % Right
    ];

    % Draw cube (camera body)
    patch('Vertices', verts', 'Faces', faces, ...
          'FaceColor', bodyColor, 'EdgeColor', 'k', 'FaceAlpha', 0.6);

    %% ===== Create Camera Lens (Cylinder) =====
    [Xc, Yc, Zc] = cylinder(0.02 * scale, 20);  % radius = 0.02
    Zc = Zc * 0.03 * scale;                     % lens depth
    lensPoints = [Xc(:)'; Yc(:)'; Zc(:)'];

    % Rotate and translate lens
    lensPoints = R * lensPoints + origin;
    Xc_new = reshape(lensPoints(1,:), size(Xc));
    Yc_new = reshape(lensPoints(2,:), size(Yc));
    Zc_new = reshape(lensPoints(3,:), size(Zc));

    % Draw lens cylinder
    surf(Xc_new, Yc_new, Zc_new, ...
         'FaceColor', lensColor, 'EdgeColor', 'none', 'FaceAlpha', 0.6);

    %% ===== Draw Front Cap Outline =====
    % Outline circle at the lens front (Z = maximum plane)
    theta = linspace(0, 2*pi, 100);
    r = 0.02 * scale;
    xc = r * cos(theta);
    yc = r * sin(theta);
    zc = ones(size(theta)) * 0.03 * scale;
    cap = R * [xc; yc; zc] + origin;

    plot3(cap(1,:), cap(2,:), cap(3,:), 'k', 'LineWidth', 0.5);
end
