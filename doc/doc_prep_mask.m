%%  Development of the yalebox-piv preprocessing masking routine.
% Masking means setting regions outside of the sand in my images to some value
% that indicates "no data", which can be ignored during PIV analysis.

%% Initialization

% setup environment
addpath('../');

% define parameters
entr_win = 11;
num_cluster = 3;

% load and original images
rgb_train = imread('doc_prep_mask_train.png');
rgb_apply = imread('doc_prep_mask_apply.png');

%% Get a manual mask for regions that are always outside the sand. 
% This step is interactive, and it is done mainly to deal with the glass support
% at the s-point. I typically also mask out obvious glare, since it is easy to
% do so, and can be helpful to segmentation down the line.

mask_manual = prep_mask_manual(rgb_train);

% plot results

%% Get parameters for automatic masking using k-means clustering
% Automatic segmentation is advantageous in that it removes the need to
% guess-and-check and makes the analysis more consistent from experiment to
% experiment.

% prepare input variables 
hsv_train = rgb2hsv(rgb_train);
hue_train = hsv_train(:,:,1);
value_train = hsv_train(:,:,3);
entropy_train = entropyfilt(value_train, true(entr_win));
 
% fit the k-means cluster model
[nr, nc] = size(hue_train);
cluster_data_train = [hue_train(:), value_train(:), entropy_train(:)];
[cluster_label, cluster_center] = kmeans(cluster_data_train, num_cluster, ...
    'Display', 'iter', ...
    'Distance', 'sqeuclidean', ...
    'Replicates', 3);
cluster_label_train = reshape(cluster_label_train, [nr, nc]);

% select the correct cluster label for the sand
figure;
imagesc(cluster_label_train);
colorbar;
while 1
    sand_label_cell = inputdlg('Sand cluster label');
    try
        sand_label = str2double(sand_label_cell{1});
    end
    if round(sand_label) == sand_label && sand_label >= 1 && sand_label <= num_cluster
        break
    end
end
close(gcf);

% shift sand cluster to label == 1
cluster_center = circshift(cluster_center, -(sand_label-1), 1); 

% plot results

%% Apply manual and automatic masking to another image from the same series

% get input layers
hsv_apply = rgb2hsv(rgb_apply);
hue_apply = hsv_apply(:,:,1);
value_apply = hsv_apply(:,:,3);
entropy_apply = entropyfilt(value_apply, true(entr_win));

% compute cluster labels for apply data
cluster_data_apply = [hue_apply(:), value_apply(:), entropy_apply(:)];
warning('off', 'stats:kmeans:FailedToConverge');
cluster_label_apply = kmeans(cluster_data_apply, num_cluster, ...
    'MaxIter', 1,...
    'Start', cluster_center);
warning( 'on', 'stats:kmeans:FailedToConverge');
cluster_label_apply = reshape(cluster_label_apply, [nr, nc]);

% generate raw mask
mask_apply = (cluster_label_apply == 1) & (mask_manual);

% plot results

%% Cleanup
% A few additional steps are needed to cleanup the data:
% * hole filling routine to take care of any anomolous holes in the sand
% * delete all but largest object (the sand)
% * remove unwanted "halo" from the upper boundary using a morphological filter
% * reapply the manual mask

% fill holes, wall off left, right and bottom
wall_lr = true(size(mask_apply, 1), 1);
mask_apply = [wall_lr, mask_apply, wall_lr];
wall_b = true(1, size(mask_apply,2));
mask_apply = [mask_apply; wall_b];
mask_apply = imfill(mask_apply, 'holes');
mask_apply = mask_apply(1:end-1, 2:end-1);

% extract largest connected object
object_label = bwlabel(mask_apply);
largest_object = mode(object_label(object_label>0));
mask_apply = object_label == largest_object;

% remove upper "halo"
disk = strel('disk', ceil(entr_win/2)+1);
mask_apply = imerode(mask_apply, disk);

% reapply manual mask
mask_apply = mask_apply & mask_manual;

% plot results

%% Apply mask to image

% sand image with red for mask
rgb_test = rgb_apply;
rgb_test(repmat(~mask_apply, [1, 1, 3])) = 0;
rgb_test(cat(3, ~mask_apply, false(size(rgb_test,1), size(rgb_test,2), 2))) = inf;

% plot results

%% Close down script

rmpath('../');