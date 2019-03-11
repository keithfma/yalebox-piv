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
labels = predict(model, features);
if iscell(labels)
    % some models (i.e., random forest) return silly cell arrays, convert to numeric
    labels = cellfun(@str2double, labels);
end
assert(all(unique(labels(:)) == [1; 2]), 'Unexpected class(es) in output labels')

% convert labels to mask
mask = labels(segments);

% TODO: use object sizes to clean things up (drop below min object size)
% TODO: use morphological filters to clean things up

% display
if ~isempty(rgb)
    hf = figure;
    hf.Name = 'Masking Results';
    
    hax1 = subplot(2, 1, 1);
    imagesc(rgb, 'AlphaData', mask == 1);
    axis equal tight
    hax1.Color = [1.0, 0.25, 0.25];
    hax1.YDir = 'normal';
    title('Sand');
    
    hax2 = subplot(2, 1, 2);
    imagesc(rgb, 'AlphaData', mask == 2);
    axis equal tight
    hax2.Color = [1.0, 0.25, 0.25];
    hax2.YDir = 'normal';
    title('Other');
    
    linkaxes([hax1, hax2], 'xy');
end