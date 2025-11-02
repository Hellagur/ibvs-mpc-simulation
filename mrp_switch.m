function ss = mrp_switch(s)
%MRP_SWICTH     Maps a Modified Rodrigues Parameter (MRP) vector to its
%               shadow set if its norm exceeds 1 to ensure numerical stability.
%
% Inputs:
%   s  - 3×1 MRP vector
%
% Outputs:
%   ss - 3×1 MRP vector, mapped to the shadow set if |s| > 1
%
% Description:
%   MRPs provide a minimal 3-parameter representation of rotation.
%   The norm of an MRP vector s satisfies |s| <= 1 in the standard set.
%   If |s| > 1, the vector is converted to its shadow set to avoid
%   singularities and maintain stable attitude computations:
%
%       s_shadow = - s / |s|


    % Compute norm of MRP vector
    s_norm = norm(s);

    % Check if shadow set mapping is needed
    if s_norm > 1
        ss = -s / s_norm;  % map to shadow set
    else
        ss = s;             % retain original
    end   
end
