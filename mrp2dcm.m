function R = mrp2dcm(s)
%MRP2DCM    Converts Modified Rodrigues Parameters (MRPs) to a Direction Cosine Matrix (DCM)
%
% Inputs:
%   s - 3×1 vector of Modified Rodrigues Parameters (MRPs)
%
% Outputs:
%   R - 3×3 rotation matrix (DCM) representing rotation from frame A to frame B
%
% Description:
%   This function computes the DCM corresponding to a given MRP vector.
%   MRPs are a minimal 3-parameter attitude representation related to
%   the rotation vector. The transformation formula is:
%
%       R = I + (8*[s]_x^2 - 4*(1 - s^T s)*[s]_x) / (1 + s^T s)^2
%
%   where [s]_x is the skew-symmetric matrix of s:
%       [s]_x = [  0   -s3   s2;
%                 s3    0   -s1;
%                -s2   s1    0 ]
%
%   The resulting DCM is normalized by its determinant to ensure it remains
%   orthonormal, improving numerical stability.


    % Construct skew-symmetric matrix from MRPs
    rodx = [   0,  -s(3),  s(2);
             s(3),     0, -s(1);
            -s(2),  s(1),    0];

    % Compute DCM from MRPs
    R = eye(3) + (8*rodx^2 - 4*(1 - s'*s)*rodx) / (1 + s'*s)^2;

    % Normalize to ensure orthonormality
    R = R / det(R);
end
