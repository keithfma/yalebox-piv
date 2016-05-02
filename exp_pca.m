% experimenting with using PCA to transform the 3-band color image to a single
% band grayscale image. The alternative is to use the hsv transform as I do now.

% conclusions: value and pca are nearly equivalent, except that pca yields a
% smoother histogram. Lab fucks everything up, it is not an option.

% considering the minor difference and the added difficulty in explaining the
% method, I plan to stick with the HSV transform.

%% local vars

mask = mask_auto & mask_manual;
rgb = rgb_train;
rgb(repmat(~mask, [1 1 3])) = 0;

%% original conversion

hsv = rgb2hsv(rgb);
hsv(repmat(~mask, [1 1 3])) = 0;
orig = hsv(:,:,3);

%% pca conversion

r = double(rgb(:,:,1));
g = double(rgb(:,:,2));
b = double(rgb(:,:,3));
data = [r(mask), g(mask), b(mask)];
[coeff,score,latent,tsquared,explained,mu] = pca(data);

new = [r(:), g(:), b(:)]*coeff(:,1);
new = reshape(new, size(rgb,1), size(rgb,2));

%% lab conversion

lab = rgb2lab(rgb);
lab(repmat(~mask, [1 1 3])) = 0;
new2 = lab(:,:,3);

%% normalize and compare

orig_norm = orig-min(orig(:));
orig_norm = orig_norm./max(orig_norm(:));

new_norm = new-min(new(:));
new_norm = new_norm./max(new_norm(:));

new2_norm = new2-min(new2(:));
new2_norm = new2_norm./max(new2_norm(:));


figure

subplot(3,1,1)
imagesc(orig_norm)
colorbar
axis equal

subplot(3,1,2)
imagesc(new_norm)
colorbar
axis equal

subplot(3,1,3)
imagesc(new2_norm)
colorbar
axis equal

linkaxes

figure
nbins = 100;

subplot(1,3,1)
hist(orig_norm(mask), nbins);

subplot(1,3,2)
hist(new_norm(mask), nbins);

subplot(1,3,3)
hist(new2_norm(mask), nbins);