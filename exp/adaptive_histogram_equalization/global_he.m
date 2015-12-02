function [eql] = global_he(im, ignore)
%
% Remap intensities of the input image so that the intensity PDF is (nearly)
% uniform.
%
% im = 2D matrix, double. Continuous-valued image matrix.
%
% ignore = Scalar, double. Values to ignore in equalization routine (masked
%   pixels).
%
% eql = 2D matrix, double. Equalized image matrix, with a uniform distribution
%   in the range [0,1]

% check inputs
validateattributes(im, {'double'}, {'2d', 'real'}, mfilename, 'im');
validateattributes(ignore, {'double'}, {'scalar', 'real'}, mfilename, 'ignore');

% identify pixels that have data
kk = im~= ignore;

% compute transform (cumulative distribution function), ignoring masked pixels
[tform_eql, tform_im] = ecdf(im(kk));

% drop repeated pixels
tform_eql = tform_eql(2:end);
tform_im = tform_im(2:end);

% apply transform (interpolate each pixel)
eql = zeros(size(im));
eql(kk) = interp1(tform_im, tform_eql, im(kk), 'linear');

% done