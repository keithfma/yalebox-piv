function mask = prep_mask_apply(model, features, segments, rgb)
% function mask = prep_mask_apply(model, features, segments, rgb)
% 
% Apply trained classifier and morphological clean-up to create a
% sand/other mask array
% 
% Arguments:
% 
%   model: ML model class, trained classifier
% 
%   features: TODO
%   
%   segments: TODO
% 
%   rgb: 3D matrix, RGB 24-bit image, optional argument, if provided
%       display resulting masked image for inspection
% 
% Returns:
%   mask: 2D matrix, 1 where sand, 0 where other
% % 

% TODO: add log messages

% set defaults
narginchk(3, 4);
if nargin < 4; rgb = []; end

% sanity checks
% TODO

% predict
fprintf('%s: predict image segment class (sand, other)\n', mfilename);
labels = predict(model, features);
if iscell(labels)
    % some models (i.e., random forest) return silly cell arrays, convert to numeric
    labels = cellfun(@str2double, labels);
end
assert(all(unique(labels(:)) == [1; 2]), 'Unexpected class(es) in output labels')

% convert labels to mask
mask = labels(segments);

% screen out by connected object area, convert to logical mask
area_threshold = 50*50;
fprintf('%s: screen out small objects (< %i pixels)\n', mfilename, area_threshold);
objects = bwlabel(mask==1) + 1;  % convert to 1-based index
object_areas = splitapply(@sum, ones(size(objects(:))), objects(:));
keep_objects = object_areas > area_threshold;
keep_objects(1) = false;  % 0 means not in original mask, drop it
mask = keep_objects(objects); 

% fill holes
fprintf('%s: fill holes\n', mfilename);
mask = imfill(mask, 'holes');

% clean up edges with morphological filters
fprintf('%s: clean up edges with morpological filters\n', mfilename);
mask = imopen(mask, strel('disk', 15));
% mask = imclose(mask, strel('disk', 15));

% display
if ~isempty(rgb)
    hf = figure;
    hf.Name = 'Masking Results';
    
    hax1 = subplot(2, 1, 1);
    imagesc(rgb, 'AlphaData', mask);
    axis equal tight
    hax1.Color = [1.0, 0.25, 0.25];
    hax1.YDir = 'normal';
    title('Sand');
    
    hax2 = subplot(2, 1, 2);
    imagesc(rgb, 'AlphaData', ~mask);
    axis equal tight
    hax2.Color = [1.0, 0.25, 0.25];
    hax2.YDir = 'normal';
    title('Other');
    
    linkaxes([hax1, hax2], 'xy');
end