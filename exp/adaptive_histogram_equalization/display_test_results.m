function [] = display_test_results(im, ignore, eql, name)
% function [] = display_test_results(im, ignore, eql, name)
%
% Display some basic results from an image equalization test.
%
% Arguments:
%
% im = 2D matrix, double. Continuous-valued image matrix.
%
% ignore = Scalar, double. Values to ignore in equalization routine (masked
%   pixels).
%
% eql = 2D matrix, double. Equalized image matrix, with a uniform distribution
%   in the range [0,1]
%
% name = String, name to include in the figure top bar to differentiate tests
%
% %

% check inputs
validateattributes(im, {'double'}, {'2d', 'real'}, mfilename, 'im');
validateattributes(ignore, {'double'}, {'scalar', 'real'}, mfilename, 'ignore');
validateattributes(eql, {'double'}, {'2d', 'real'}, mfilename, 'im');
validateattributes(name, {'char'}, {}, mfilename, 'name');

% plot original and equalized image
figure
set(gcf, 'Name', name);

subplot(2,1,1)
imagesc(im);
caxis([0,1]);
axis off
title('Original grayscale');

subplot(2,1,2)
imagesc(eql);
caxis([0,1]);
axis off
title('Equalized');

% plot empirical CDF and PDF for original and equalized images
figure
set(gcf, 'Name', name);

roi = im~=0;

[cdf, val] = ecdf(im(roi));
pdf = diff(cdf)./diff(val);

subplot(2,2,1)
plot(val, cdf, 'Marker', '.');
title('Original CDF');

subplot(2,2,3)
plot(val(2:end), pdf, 'Marker', '.');
title('Original PDF');

[cdf, val] = ecdf(eql(roi));
pdf = diff(cdf)./diff(val);

subplot(2,2,2)
plot(val, cdf, 'Marker', '.');
title('Equalized CDF');

subplot(2,2,4)
plot(val(2:end), pdf, 'Marker', '.');
title('Equalized PDF');