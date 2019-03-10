function segments = prep_segment(img, scale, sigma, min_size)
% function segments = prep_segment(img, scale, sigma, min_size)
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
%   img: TODO
%   scale: TODO
%   sigma: TODO
%   min_size: TODO
%
% Returns:
%   TODO
% %

% set defaults
narginchk(1, 4);
if nargin < 2; scale = 200; end
if nargin < 3; sigma = 0.5; end
if nargin < 4; min_size = 50; end

% get path to executable
[~, python_interpreter, ~] = pyversion();
python_script = [mfilename('fullpath'), '.py'];

% create temporary files for IO
input_file = [tempname, '.mat'];
output_file = [tempname, '.mat'];
cleaner = onCleanup(@() delete(input_file, output_file));

% save image to mat file
save(input_file, 'img');

% run python script
cmd = sprintf('%s %s %s %s --scale %i --sigma %.03f --min_size %i', ...
    python_interpreter, python_script, input_file, output_file, scale, sigma, min_size);
status = system(cmd);

% unpack results
results = load(output_file, 'segments');
segments = results.segments;