function features = prep_mask_features(rgb)
% function features = prep_mask_features(rgb)
% 
% Create feature matrix (rows as observations, columns as features) from
% the input experiment image
% 
% Arguments:
%   rgb: 3D matrix, RGB 24-bit image
% 
% Returns:
%   features: 2D matrix, rows as observations (pixels)
% % 
hsv = rgb2hsv(rgb);

entropy_kernel = strel('disk', 5).Neighborhood;
h_e = entropyfilt(hsv(:, :, 1), entropy_kernel);
s_e = entropyfilt(hsv(:, :, 2), entropy_kernel);
v_e = entropyfilt(hsv(:, :, 3), entropy_kernel);

as_feature = @(x) double(reshape(x, [], 1));

features = [...
    as_feature(hsv(:, :, 1)), ...
    as_feature(hsv(:, :, 2)), ...
    as_feature(hsv(:, :, 3)), ...
    as_feature(h_e), ...
    as_feature(s_e), ...
    as_feature(v_e), ...
    ];