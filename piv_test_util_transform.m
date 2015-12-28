function [x_out, y_out] = piv_test_util_transform(tform, x_in, y_in, fwd)
% function [x_out, y_out] = piv_test_util_transform(tform, x_in, y_in, fwd)
% 
% tform = 2x3 affine transformation matrix
%
% x_in, y_in = vectors of x, y point coordinates
%
% fwd = Scalar, logical, flag indicating if the transform should be forward (1)
%   or reverse (0)
%
% x_out, y_out = vectors of transformed x,y point coordinates
%
% %

% get full transform matrix
A = [tform; 0 0 1]; 
if ~fwd; 
    A = inv(A); 
end

% transform points, maintaining vector shape
pts_in = [x_in(:)'; y_in(:)'; ones(1, length(x_in))];

pts_out = A*pts_in;

x_out = reshape(pts_out(1,:), size(x_in));
y_out = reshape(pts_out(2,:), size(y_in));

end