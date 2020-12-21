% Modify the prepped image file Mark created to match yalebox-piv v0.7.1 expected format
% %

% note: assumes you ran update_path('main') in the root of yalebox_piv, so that the update_path
%   function is available to run here
update_path('newmatic');

SRC_FILE_NAME = 'images_src.mat';
REF_FILE_NAME = 'images_ref.mat';
DEST_FILE_NAME = 'images.mat';

% read coordinates from source and get array sizes
% note: we use the destination variable names everywhere except when accessing the source file
% note: the RGB images in the src file are also padded, which means we must use the "extended"
%   coordinates for all our arrays (i.e., x == x_ext)
% note: src file has units of mm, but we use meters and so convert
src_mat = matfile(SRC_FILE_NAME, 'Writable', false);
ref_mat = matfile(REF_FILE_NAME, 'Writable', false);
MM_TO_M = 1/1000;

x_ext = src_mat.xPad*MM_TO_M;
num_x_ext = length(x_ext);

y_ext = src_mat.yPad*MM_TO_M;
num_y_ext = length(y_ext);

band = ref_mat.band;  % not present in src file
num_band = length(band);

num_img =  size(src_mat, 'imgGray', 3);


% declare the output file, properly chunked
delete(DEST_FILE_NAME);
dest_mat = newmatic(...
    DEST_FILE_NAME, ...
    newmatic_variable('img', 'uint8', [num_y_ext, num_x_ext, num_band, num_img], [num_y_ext, num_x_ext, num_band, 1]), ...
    newmatic_variable('img_ext', 'single', [num_y_ext, num_x_ext, num_img], [num_y_ext, num_x_ext, 1]), ...
    newmatic_variable('mask', 'logical', [num_y_ext, num_x_ext, num_img], [num_y_ext, num_x_ext, 1]), ...
    newmatic_variable('mask_ext', 'logical', [num_y_ext, num_x_ext, num_img], [num_y_ext, num_x_ext, 1]) ...
);

% populate the output file one variable at a time
fprintf('write coordinates\n');
dest_mat.band = band;
dest_mat.x = x_ext;  % using _ext is right
dest_mat.x_ext = x_ext;
dest_mat.y = y_ext; % using _ext is right
dest_mat.y_ext = y_ext;
dest_mat.step = 0:(num_img-1);

fprintf('write RGB images\n');
% note: must write to slices, or else chunking is clobbered
tmp = src_mat.img;
for ii = 1:size(tmp, 4)
    dest_mat.img(:, :, :, ii) = tmp(:, :, :, ii);
end

fprintf('write grayscale images\n');
% note: must write to slices, or else chunking is clobbered
% note: must cast uint8 to single
tmp = single(src_mat.imgGray)/255;
for ii = 1:size(tmp, 3)
    dest_mat.img_ext(:, :, ii) = tmp(:, :, ii);
end

fprintf('write masks\n');
% note: must write to slices, or else chunking is clobbered
tmp = src_mat.mask;
for ii = 1:size(tmp, 3)
    dest_mat.mask(:, :, ii) = tmp(:, :, ii);
    dest_mat.mask_ext(:, :, ii) = tmp(:, :, ii);  % only one mask in src
end
