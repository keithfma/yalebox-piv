function [color_ref, subset] = get_color_ref(image_file, subset)
% function [color_ref, subset] = get_color_ref(image_file, subset)
%
%
% Return the mean color (24-bit RGB) of a subset of the the image in
% image_file. The subset region can be selected interactively (only
% image_file is provided) or can be speficied in the function call (both
% image_file and subset are provided).
%
% Arguments:
%
%   image_file = String. Filename of the image to be analyzed.
%   subset = 4-element vector. Limits of the subset to be averaged, in the
%           form [xmin ymin width height]. Will be set interactively if not 
%           provided as an input argument.
%   color_ref = 3 element vector. 24-bit RGB representation of the average
%               color of sand.
%
% Keith Ma, July 2015

%% read in data, check for sane inputs

assert(exist(image_file, 'file') == 2, 'image_file is not a valid filename');
img = imread(image_file);
assert(isa(img, 'uint8') && size(img,3) == 3, 'image is not 24-bit RGB');

if nargin == 2
    tmp = size(subset);
    assert(numel(tmp) == 2 && min(tmp) == 1 && max(tmp) == 4, ...
        'subset is not a 4-element vector');
    assert(subset(1) >= 1 && subset(1)+subset(3) <= size(img,2) ...
        && subset(2) >= 1 && subset(2)+subset(4) <= size(img,1), ...
        'subset extends beyond the limits of the image data');
end

%% get subset extent (if not supplied)

if nargin ==1
    while 1
        [img_subset, subset] = imcrop(img);
        hfig = figure;
        imshow(img_subset);
        button = questdlg('Continue using this region?','get_color_ref','Continue','Retry','Cancel', 'Retry');
        close(hfig);
        switch button
            case 'Continue'                
                break
            case 'Retry'
                continue
            case 'Cancel'
                error('get_color_ref canceled by user');
        end
    end
else
    [img_subset, subset] = imcrop(img, subset);
end

%% compute average value

color_ref = [0, 0, 0];
for i = 1:3
    color_ref(i) = mean(mean(double(img_subset(:,:,i))));
end