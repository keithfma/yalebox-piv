% Script. Explore imposing lower boundary condition by padding the initial and
% final  images.
%

clear

%% v == 0
% The pad varies (randomly) in the y-direction, constant in the x-direction, the
% same in both images. Adding this pad should sharpen the peak in the
% y-direction


% load test data (samp, intr, xcr):
load simple_shear_upper_bnd.mat

% add pad 
pad_val = rand(max(size(intr)), 1)*range(intr(:));

samp_pad0 = find(all(samp == 0, 2), 1, 'first');
samp_npad = size(samp, 1)-samp_pad0+1;
samp_pad = repmat(pad_val(1:samp_npad), 1, size(samp, 2));
samp(samp_pad0:end, :) = samp_pad;

intr_pad0 = find(all(intr == 0, 2), 1, 'first');
intr_npad = size(intr, 1)-intr_pad0+1;
intr_pad = repmat(pad_val(1:intr_npad), 1, size(intr, 2));
intr(intr_pad0:end, :) = intr_pad;

% redo cross-correlation
pad_xcr = normxcorr2(samp, intr);

%% du/dy == 0 outside the image domain
% The pad copies out the last row in each image

% load test data (samp, intr, xcr):
load simple_shear_upper_bnd.mat

samp_pad0 = find(all(samp == 0, 2), 1, 'first');
samp_npad = size(samp, 1)-samp_pad0+1;
samp_pad = repmat(samp(samp_pad0-1, :), samp_npad, 1);
samp(samp_pad0:end, :) = samp_pad;

intr_pad0 = find(all(intr == 0, 2), 1, 'first');
intr_npad = size(intr, 1)-intr_pad0+1;
intr_pad = repmat(intr(intr_pad0-1, :), intr_npad, 1);
intr(intr_pad0:end, :) = intr_pad;

% redo cross-correlation
pad_xcr = normxcorr2(samp, intr);

%% both
% The pad copies out the last row, and adds a random pattern in the y-direction

% load test data (samp, intr, xcr):
load simple_shear_upper_bnd.mat

% add pad 
pad_val = rand(max(size(intr)), 1)*range(intr(:));
pad_val = pad_val - mean(pad_val(:));

samp_pad0 = find(all(samp == 0, 2), 1, 'first');
samp_npad = size(samp, 1)-samp_pad0+1;
samp_pad_a = repmat(samp(samp_pad0-1, :), samp_npad, 1);
samp_pad_b = repmat(pad_val(1:samp_npad), 1, size(samp, 2));
samp(samp_pad0:end, :) = samp_pad_a+samp_pad_b;

intr_pad0 = find(all(intr == 0, 2), 1, 'first');
intr_npad = size(intr, 1)-intr_pad0+1;
intr_pad_a = repmat(intr(intr_pad0-1, :), intr_npad, 1);
intr_pad_b = repmat(pad_val(1:intr_npad), 1, size(intr, 2));
intr(intr_pad0:end, :) = intr_pad_a+intr_pad_b;

% redo cross-correlation
pad_xcr = normxcorr2(samp, intr);
