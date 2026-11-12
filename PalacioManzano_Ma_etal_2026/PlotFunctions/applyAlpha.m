function rgbA = applyAlpha(rgb, alpha, bg_rgb)
% APPLYALPHA Blends an RGB triplet with a background color using an alpha value.
%
%   rgb_out = applyAlpha(rgb_in, alpha)
%   rgb_out = applyAlpha(rgb_in, alpha, bg_rgb)
%
%   Inputs:
%       rgb_in - 1x3 vector [R, G, B] (values between 0 and 1)
%       alpha  - Scalar transparency value (0 = fully transparent, 1 = fully opaque)
%       bg_rgb - (Optional) 1x3 background RGB triplet. Defaults to [1 1 1] (white).
%
%   Output:
%       rgb_out - 1x3 blended RGB triplet [R, G, B]

    if nargin < 3 || isempty(bg_rgb)
        bg_rgb = [1, 1, 1]; % Default background is white
    end

    % Input validation and bounds checking
    validateattributes(rgb, {'numeric'}, {'vector', 'numel', 3, '>=', 0, '<=', 1});
    validateattributes(alpha,  {'numeric'}, {'scalar', '>=', 0, '<=', 1});
    validateattributes(bg_rgb, {'numeric'}, {'vector', 'numel', 3, '>=', 0, '<=', 1});

    % Ensure row vectors
    rgb = rgb(:)';
    bg_rgb = bg_rgb(:)';

    % Alpha blending equation: C_out = alpha * C_src + (1 - alpha) * C_bg
    rgbA = alpha * rgb + (1 - alpha) * bg_rgb;
end