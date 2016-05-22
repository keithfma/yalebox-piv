function [fimg, fy] = movie_frame_flip(img, y)
%
% Flip the veritical axis of the image so that the image top (i.e. maximum
% y value) is garunteed to lie at row == 1.
%
% Arguments:
%
% 
%
% Keith Ma, August 2015

if y(1) < y(end)
    fy = y(end:-1:1);
    fimg = flipud(img);
else 
    fy = y;
    fimg = img;
end