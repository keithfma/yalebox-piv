function [img_fill, mask_fill] = prep_fill(img, mask, method, show)
% function img_fill = prep_fill(img, mask, method, show)
% 
% Fill non-sand regions of the input image to enable PIV on full image domain
%
% Arguments:
%   img: 2D matrix, rectified/cropped, and padded grayscale image
%   mask: 2D matrix, logical mask of sand pixels in img
%   method: String, name of method to use for filling pixels, current options are: "none", default
%       is "none"
%   show: Logical, set true to display results for visual debugging, default is false
% 
% Outputs:
%   img_fill: 2D matrix, filled grayscale image
%   mask_fill: 2D matrix, image mask matching size of img_fill
% %

% FIXME: is this note still true?
% note: x-direction is padded with zeros and left masked. There is no
%   sane way to pad regions where sand grains enter or leave the frame.

% set defaults
if nargin < 3; method = 'none'; end
if nargin < 4; show = false; end 

% sanity checks
validateattributes(img, {'double', 'single'}, {'2d'});
validateattributes(mask, {'logical'}, {'2d', 'size', size(img)});
validateattributes(method, {'char'}, {});
validateattributes(show, {'logical'}, {'scalar'});

fprintf('%s: fill image using method: "%s"\n', mfilename, method); 

if strcmp(method, 'none')
    % do nothing, and use original mask
    img_fill = img;
    mask_fill = mask;
    
elseif strcmp(method, 'mirror')
    % mirror vertically across sand boundary
    warning('fill method "%s" known to produce incorrect results', method);
    [img_fill, mask_fill] = prep_fill_mirror(img, mask);
    
else
    error('Unknown fill method "%s"', method);

end

if show
    hf = figure;
    hf.Units = 'normalized';
    hf.Position = [0, hf.Position(2), 1, hf.Position(4)];
    
    imagesc(img_fill);
    hold on
    bnds = bwboundaries(mask);  % not fill, original mask
    for ii = 1:length(bnds)
        bnd = bnds{ii};
        plot(bnd(:, 2), bnd(:, 1), 'r', 'LineWidth', 2);
    end
            
    colormap('gray');
    hax = gca;
    hax.YDir = 'normal';
    hax.YTick = [];
    hax.XTick = [];
    axis equal tight
    title(sprintf('Image filled with method "%s"', method'));
end