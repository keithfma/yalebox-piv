function segments = prep_mask_segments(img, scale, sigma, min_size, show)
% function segments = prep_mask_segments(img, scale, sigma, min_size, show)
%
% Multiband Felzenswalb image segmentation
%
% Note: implemented by calling a python module to use skimage routines,
%   this is awkward, but no MATLAB implementation exists at present. Some
%   effort is taken to make this work on different platforms, but it is
%   assumed that the active python environment includes scipy and
%   scikit-image
%
% Note: uses the python interpreter that MATLAB from the current MATLAB
%   configuration, which you can see by running `pyversion`, and can change
%   by running `pyversion path/to/python/executable`
%
% Arguments:
%   img: N-D multiband image to be segmented
%   scale: TODO
%   sigma: TODO
%   min_size: TODO
%   show: display results for smell test, assumes first 3 bands of the
%       input img with plot nicely
%
% Returns:
%   integer labels for image segments, footprint is same as first two
%   dimensions of the input img
% %

% set defaults
narginchk(1, 5);
if nargin < 2; scale = 200; end
if nargin < 3; sigma = 0.5; end
if nargin < 4; min_size = 50; end
if nargin < 5; show = false; end

% get path to executable
[~, python_interpreter, ~] = pyversion();
python_script = [mfilename('fullpath'), '.py'];

% create temporary files for IO
input_file = [tempname, '.mat'];
output_file = [tempname, '.mat'];
cleaner = onCleanup(@() delete(input_file, output_file));

% save image to mat file
fprintf('%s: write input image to disk\n', mfilename);
save(input_file, 'img');

% run python script
fprintf('%s: segment input image\n', mfilename);
cmd = sprintf('%s %s %s %s --scale %i --sigma %.03f --min_size %i', ...
    python_interpreter, python_script, input_file, output_file, scale, sigma, min_size);
status = system(cmd);
assert(status == 0, 'External command failed: %s', cmd);

% unpack results
fprintf('%s: load output image from disk\n', mfilename);
results = load(output_file, 'segments');
segments = results.segments + 1;  % shift to 1-based index, helps with downstream compatibility

% optionally display results
if show
    figure;
    disp_img = imoverlay(...
        img(:,:,3), ...
        boundarymask(double(segments)), ...
        'cyan'); 
    imshow(disp_img);
    title(sprintf(...
        'Segmentation boundaries for scale=%i, sigma=%.03f, min\_size=%i', ...
        scale, sigma, min_size));
    pause
end