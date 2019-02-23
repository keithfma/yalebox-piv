function unpacked = unpack_polygons(packed)
% function unpacked = unpack_polygons(packed)
% 
% Arguments:
%   packed = polygon vertices as 2d array with points in columns and NaNs
%       separating each polygon
% 
% Returns:
%   unpacked = cell array or 2D matrices, each cell contains vertices for
%       one polygon with x-coord in row 1 and y-coord in row 2
% % 

unpacked = {};    
ini = 1;
while ini < size(packed, 2)
    fin = find(isnan(packed(1, ini:end)), 1, 'first') + ini - 1;
    unpacked{end+1} = packed(:, ini:(fin-1))'; %#ok!
    ini = fin + 1;
    if ini >= size(packed, 2)
        break
    end
end