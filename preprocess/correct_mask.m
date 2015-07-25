function [mask] = correct_mask(img, entropy_radius, open_radius, bw_threshold)
% function [mask] = correct_mask(img, entropy_radius, open_radius, bw_threshold)
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
%   img = Double, 3D matrix. RGB image, normalized to the range [0, 1]
%   entropy_radius = Scalar. Radius of neighboorhood included in the
%       entropy filter calculation for each pixel, [pixels]
%   open_radius = Scalar. Radius of the structuring element used for the
%       morphological opening filter, [pixels]
%   bw_threshold = Scalar in the range [0,1], threshold value used to
%       convert mask image from grayscale to black-and-white.
%   mask = Logical, 2D matrix. TRUE where there is sand, FALSE elsewhere
%
% Keith Ma, July 2015

%% read in data, check for sane inputs

assert(isa(img, 'double'), 'img is not of type double');
assert(size(img,3) == 3, 'img is not a 3D array (RGB)');
assert(max(img(:)) <= 1 && min(img(:)) >= 0, 'img is not in the range [0,1]');

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
entr = entropyfilt(img, nhood);

keyboard

s = strel('disk', open_radius);
img = imopen(img, s);

img = img-min(img(:));
img = img./max(img(:));
mask = im2bw(img, bw_threshold);
