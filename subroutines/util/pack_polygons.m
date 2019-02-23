function packed = pack_polygons(unpacked)
% function packed = pack_polygons(unpacked)
% 
% Arguments:
%   unpacked = cell array or 2D matrices, each cell contains vertices for
%       one polygon with x-coord in row 1 and y-coord in row 2
% 
% Returns:
%   packed = polygon vertices as 2d array with points in columns and NaNs
%       separating each polygon
% % 

packed = [];
for ii = 1:length(unpacked)
   this_poly = unpacked{ii};
   packed = [packed, this_poly', nan(2,1)]; %#ok!
end
