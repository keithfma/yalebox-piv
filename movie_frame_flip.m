function [fimg, fy] = movie_frame_flip(img, y)
% function [fimg, fy] = movie_frame_flip(img, y)
%
% Flip the veritical axis of the image so that the image top (i.e. maximum
% y value) is guaranteed to lie at row == 1.
%
% Arguments:
%   img: Image matrix, with x-direction along columns, and y-direction along rows
%   y: Vector, y-direction coordinates corresponding to img
%   fimg: Flipped image matrix
%   fy: Flipped coordinate vector, corresponding to fimg
% 
% % Keith Ma

if y(1) < y(end)
    fy = y(end:-1:1);
    fimg = flipud(img);
else 
    fy = y;
    fimg = img;
end
