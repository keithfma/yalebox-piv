% Explore clustering algorithms for image segmentation

%% load image and generate possible input layers

% parameters
image_file = 'test.png';
entr_win = 11;
med_win = 25;

% possible input layers
rgb = imread(image_file);
rgb_red = rgb(:,:,1);
rgb_green = rgb(:,:,2);
rgb_blue = rgb(:,:,3);
hsv = rgb2hsv(rgb);
hsv_hue = hsv(:,:,1);
hsv_saturation = hsv(:,:,2);
hsv_value = hsv(:,:,3);
[col, row] = meshgrid(1:size(rgb,2), 1:size(rgb,1));
entropy = entropyfilt(hsv_value, true(entr_win));
cform = makecform('srgb2lab');
lab = applycform(rgb,cform);
lab_light = lab(:,:,1);
lab_a = lab(:,:,2);
lab_b = lab(:,:,3);
hsv_value_median = medfilt2(hsv_value, [med_win, med_win]);
hsv_hue_median = medfilt2(hsv_hue, [med_win, med_win]);

% best layers are:
% ...hsv_hue
% ...hsv_value
% ...entropy  

% ...hsv_saturation XXX
% ...lab_light XXX
% ...row XXX
% ...col XXX

%% try kmeans clustering with one variable only

var = entropy;
num_cluster = 3;

data = reshape(var, numel(var), 1);
label = kmeans(data, num_cluster); 
label = reshape(label, size(var));

% remarkably successful!

%% try kmeans clustering with multiple variables

[nr, nc] = size(rgb_red);
num_cluster = 3;

data = [hsv_hue(:), hsv_value(:), hsv_saturation(:), entropy(:)];
label = kmeans(data, num_cluster); 
label = reshape(label, [nr, nc]);
imagesc(label);

% close to done...

%% try again with kmeans clustering

[nr, nc] = size(rgb_red);
num_cluster = 5;

data = [hsv_hue_median(:), hsv_value_median(:)];
label = kmeans(data, num_cluster); 
label = reshape(label, [nr, nc]);
imagesc(label);

%% try EM clustering with one multiple variables

[nr, nc] = size(rgb_red);
num_cluster = 5;

data = [hsv_hue(:), hsv_value(:), entropy(:)];
gmfit = fitgmdist(data, num_cluster, ...
    'CovarianceType', 'diagonal', ...
    'SharedCovariance', false, ...
    'RegularizationValue',0.01);
label = cluster(gmfit, data);
label = reshape(label, [nr, nc]);
imagesc(label);

% best is above...but very slow to compute model...and overshoots a bit due to entropy