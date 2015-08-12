function [fimg, fy] = yalebox_movie_aux_flip(img, y)
%
% Flip the veritical axis of the image so that the image top (i.e. maximum
% y value) is garunteed to lie at row == 1.
%
% Arguments:
%
% Keith Ma, August 2015

% debug
fimg = img;
fy = y;