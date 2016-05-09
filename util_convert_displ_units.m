function [xx, yy, uu, vv, mm, bbox] = ...
    util_convert_displ_units(xx, yy, uu, vv, mm, coord, displ, bbox)
%
% Convert units of displacement data and its coordinates. Conversion may be
% a constant factor (e.g. m -> cm) or data dependent (e.g. normalize by
% incoming section displacement).
%
% Note on normalized displacements: displacement vector components and are
% divided by the median displacement magnitude within the specified
% bounding box. This is useful to facilitate comparison between steps of
% unequal size. Inter-step inconsistency is a shortcoming of the
% experimental apparatus. The normalization factor is a user-selected
% characteristic velocity, computed as the absolute value of the median of
% the velocity mangitude in a user-defined window. For most sandbox
% experiments, the logical choice is the incoming section.
%
% Arguments:
%
% xx, yy = 
%
% uu, vv, mm = 
%
% coord = 
%
% displ =
%
% bbox = 
% %

% convert coordinate units
switch coord
    
    case 'm'
        % do nothing        
        
    case 'cm'
        xx = xx*100;
        yy = yy*100;
        
end

% convert displacement units
switch displ
    
    case 'm/step'
        % do nothing
    
    case 'mm/step'
        uu = uu*1000;
        vv = vv*1000;
        mm = mm*1000;
        
    case '1'
        if nargin < 8 || isempty(bbox)
            figure;
            imagesc(xx, yy, mm, 'AlphaData', ~isnan(mm));
            set(gca, 'YDir', 'Normal');
            title('Select bounding box for normalization using mouse');
            bbox = getrect();
            close(gcf);
        end
        
        % get normalization factor
        x_in_bbox = (xx >= bbox(1)) & (xx <= (bbox(1)+bbox(3)));
        y_in_bbox = (yy >= bbox(2)) & (yy <= (bbox(2)+bbox(4)));
        in_bbox = logical(bsxfun(@times, x_in_bbox(:)', y_in_bbox(:)));
        norm = nanmedian(mm(in_bbox));
        
        % apply normalization
        uu = uu./norm;
        vv = vv./norm;
        mm = mm./norm;
               
end