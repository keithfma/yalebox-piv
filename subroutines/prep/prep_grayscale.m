function gray = prep_grayscale(rgb)
% function gray = prep_grayscale(rgb)
%
% Covert 24-bit RGB to grayscale (float with range 0-1)
%
% Arguments
%   rgb: 24-bit RGB image
%   gray: grayscale version of "rgb"
% %

% TODO: explore alternative methods to generate grayscale that maximize information content
%   of the resulting image (perceptual accuracy is not so important here)

validateattributes(rgb, {'uint8'}, {'ndims', 3});

hsv = rgb2hsv(rgb);
gray = hsv(:, :, 3);
