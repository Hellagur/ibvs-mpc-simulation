function Plot_arrow3D(start_pt, end_pt, color, alpha, width, head_length)
%PLOT_ARROW3D   Draw a 3D directional arrow between two points.
%
% This function draws a smooth 3D arrow (shaft + cone head) between 
% the specified start and end points, which can be used to visualize 
% directions, forces, velocities, or camera frames in space.
%
% Inputs:
%   start_pt     — 1×3 vector, arrow start point [x, y, z].
%   end_pt       — 1×3 vector, arrow end point [x, y, z].
%   color        — 1×3 RGB vector specifying arrow color, e.g., [0, 0, 1].
%   alpha        — scalar in [0,1], transparency level (1 = opaque).
%   width        — (optional) shaft radius (default: 0.05).
%   head_length  — (optional) cone head length (default: 0.2).
%
% Example:
%   Plot_arrow3D([0,0,0], [1,1,1], [1,0,0], 0.8, 0.03, 0.15);


    % Default parameters
    if nargin < 5, width = 0.05; end
    if nargin < 6, head_length = 0.2; end

    % Compute direction vector and length
    dir = end_pt - start_pt;
    L = norm(dir);
    if L == 0, return; end
    dir = dir / L;

    % Define shaft length (excluding the cone head)
    shaft_len = L - head_length;
    
    %% ========== Draw Shaft (Cylinder) ==========
    % Create cylinder mesh (radius = width, height = shaft_len)
    [xc, yc, zc] = cylinder(width, 12);
    zc = zc * shaft_len;

    % Convert surface mesh to patch format for 3D transformation
    shaft = surf2patch(xc, yc, zc);

    % Rotate and translate the shaft to align with direction vector
    shaft_verts = transform_arrow(shaft.vertices, dir, start_pt);

    % Render the shaft
    patch('Faces', shaft.faces, 'Vertices', shaft_verts, ...
          'FaceColor', color, 'EdgeColor', 'none', 'FaceAlpha', alpha);

    %% ========== Draw Head (Cone) ==========
    % Create cone mesh (base radius = 2×shaft width)
    [xh, yh, zh] = cylinder([width*2.0, 0], 12);
    zh = zh * head_length;

    % Convert cone surface to patch format
    cone = surf2patch(xh, yh, zh);

    % Rotate and translate the cone to the end of the shaft
    cone_verts = transform_arrow(cone.vertices, dir, start_pt + shaft_len * dir);

    % Render the cone head
    patch('Faces', cone.faces, 'Vertices', cone_verts, ...
          'FaceColor', color, 'EdgeColor', 'none', 'FaceAlpha', alpha);
end


function verts_out = transform_arrow(verts, direction, origin)
%TRANSFORM_ARROW Rotate and translate arrow vertices to align with a given direction.
%
% This helper function aligns the 3D geometry of the arrow 
% (originally along the z-axis) with an arbitrary direction vector 
% and translates it to the specified origin.
%
% Inputs:
%   verts     — Nx3 vertex matrix of the object (in local coordinates).
%   direction — 3×1 unit vector specifying desired arrow direction.
%   origin    — 1×3 vector, new origin of the arrow base.
%
% Output:
%   verts_out — Nx3 transformed vertices in world coordinates.

    % Normalize the direction vector
    z = [0 0 1]';
    dir = direction(:) / norm(direction);

    % Compute rotation axis and angle using Rodrigues' formula
    v = cross(z, dir);    % rotation axis
    s = norm(v);          % sine of rotation angle
    c = dot(z, dir);      % cosine of rotation angle

    % Construct rotation matrix
    if s == 0
        % Case: direction is parallel or anti-parallel to z-axis
        R = eye(3) * sign(c);
    else
        % Skew-symmetric cross-product matrix
        vx = [  0   -v(3)  v(2);
               v(3)   0   -v(1);
              -v(2)  v(1)   0 ];
        % Rodrigues rotation formula
        R = eye(3) + vx + vx^2 * ((1 - c) / (s^2));
    end

    % Apply rotation and translation
    verts_out = (R * verts')' + origin;
end
