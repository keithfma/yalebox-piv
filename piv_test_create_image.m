function [ini, fin, ini_roi, fin_roi, xx, yy] = ...
    piv_test_create_image(tform, bnd_mean, bnd_ampl, bnd_freq)
%
% Create a synthetic image pair by deforming a template image with a homogenous
% deformation + constant offset, and imposing a sinusoidal boundary on the
% initial image. The image size is chosen so that both the initial and final
% (deformed) images are fully populated, and thus depends on both the
% transformation and the size of the image template.

%% initialize

% set defaults
narginchk(0,4);
if nargin == 0 || isempty(tform)  
    tform = [1, 0.05,  5; 0,   1,  5]; 
end 
if nargin < 2 || isempty(bnd_mean)
    bnd_mean = 0.75;
end
if nargin < 3 || isempty(bnd_ampl)
    bnd_ampl = 0.1;
end
if nargin < 4 || isempty(bnd_freq) 
    bnd_freq = 1;
end

% check for sane inputs
validateattributes(tform, {'numeric'}, {'2d', 'size', [2, 3]});
validateattributes(bnd_mean, {'numeric'}, {'scalar'});
validateattributes(bnd_ampl, {'numeric'}, {'scalar'});
validateattributes(bnd_freq, {'numeric'}, {'scalar'});

