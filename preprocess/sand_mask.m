function [mask] = sand_mask(img, entropy_radius, open_radius, bw_threshold)
% function [mask] = sand_mask(img, entropy_radius, open_radius, bw_threshold)
%
% Create a logical mask for a color image (image_file) that is TRUE where
% there is sand and FALSE elsewhere. This can be used to remove (set to 0)
% the background in a image prior to PIV analysis or other applications.
%
% Sand is identified by filtering and thresholding the image data. The
% steps are:
% 
%   (1) Convert color to grayscale
%   (2) Apply entropy filter, which provides a measure of local texture
%   (3) Apply opening filter, which eliminates narrow curvilinear objects
%       typical of reflections
%   (4) Convert to a T/F mask by simple thresholding
%
% Arguments:
%   img = MxNx3 matrix of type uint8, a 24-bit (Truecolor) RGB image 
%   entropy_radius = Scalar. Radius of neighboorhood included in the
%       entropy filter calculation for each pixel, [pixels]
%   open_radius = Scalar. Radius of the structuring element used for the
%       morphological opening filter, [pixels]
%   bw_threshold = Scalar in the range [0,1], threshold value used to
%       convert mask image from grayscale to black-and-white.
%
% Keith Ma, July 2015

% % DEBUG: Run as a script with predefined inputs
% img = imread('test_data/test_a_00_crop.png');
% entropy_radius = 9; % pixels
% open_radius = 36; % pixels
% bw_threshold = 0.75; % normalized

%% read in data, check for sane inputs

assert(isa(img, 'uint8') && size(img,3) == 3, 'image is not 24-bit RGB');

assert(numel(entropy_radius) == 1, 'entropy_radius is not a scalar');
assert(entropy_radius >= 0, 'entropy_radius is not positive');

assert(numel(open_radius) == 1, 'open_radius is not a scalar');
assert(open_radius >= 0, 'open_radius is not positive');

assert(numel(bw_threshold) == 1, 'bw_threshold is not a scalar');
assert(bw_threshold >= 0 && bw_threshold <= 1, ...
    'bw_threshold is not in the range [0, 1]');


%% compute mask

img = rgb2gray(img);

nhood = getnhood(strel('disk', entropy_radius)); 
img = entropyfilt(img, nhood);

s = strel('disk', open_radius);
img = imopen(img, s);

img = img-min(img(:));
img = img./max(img(:));
mask = im2bw(img, bw_threshold);
