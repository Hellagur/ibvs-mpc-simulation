function s = dcm2mrp(R)
%DCM2MRP    Converts a Direction Cosine Matrix (DCM) to Modified Rodrigues Parameters (MRPs)
%
% Inputs:
%   R - 3×3 rotation matrix (DCM) representing rotation from frame A to frame B
%
% Outputs:
%   s - 3×1 vector of Modified Rodrigues Parameters (MRPs)
%
% Description:
%   This function computes the MRPs corresponding to a given rotation matrix.
%   The MRPs provide a minimal 3-parameter representation of attitude,
%   related to the classical rotation vector by a stereographic projection.
%
%   The formula used is:
%       s = (1 / (sqrt(trace(R)+1) * (sqrt(trace(R)+1)+2))) * ...
%           [R(2,3)-R(3,2); R(3,1)-R(1,3); R(1,2)-R(2,1)]
%
%   If the resulting MRP norm is greater than 1, a shadow set transformation
%   is applied to ensure |s| <= 1. This improves numerical stability.
%

    % Compute intermediate scalar
    x = sqrt(trace(R) + 1);

    % Compute MRPs from DCM
    s = [R(2,3) - R(3,2);
         R(3,1) - R(1,3);
         R(1,2) - R(2,1)] / (x * (x + 2));

    % Apply shadow set if norm exceeds 1
    s = mrp_switch(s);
end
