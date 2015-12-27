function [ini, fin, ini_roi, fin_roi, xx, yy] = piv_test_create_image()
%
% Create a synthetic image pair by deforming a template image with a homogenous
% deformation + constant offset, and imposing a sinusoidal boundary on the
% initial image.
%