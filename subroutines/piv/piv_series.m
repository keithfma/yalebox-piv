function [] = piv_series(...
    output_file, input_file, step_range, gap, samp_len, samp_spc, ...
    intr_len, num_pass, valid_radius, valid_max, valid_eps, ...
    min_frac_data, min_frac_overlap, notes)
% function [] = piv_series(...
%     output_file, input_file, step_range, gap, samp_len, samp_spc, ...
%     intr_len, num_pass, valid_radius, valid_max, valid_eps, ...
%     min_frac_data, min_frac_overlap, notes)
%
% Run PIV analysis for a given input series. Input image data  expected to
% be a MAT file as created by prep_series(). Results are saved in a new
% MAT file which includes all relevant metadata.
%
% Arguments:
%
%   output_file: String, name of the output MAT file containing PIV results
%   input_file: String, name of the input MAT file containing pre-processed
%       image data
%   step_range: 2-element vector, [initial, final] steps of input images series
%       to include in PIV analysis. Either element can be set to NaN to use the
%       full range.
%   gap : Scalar, integer, gap between analyzed images in steps, for
%       example, for an initial image at step 3, setting gap -> 1 would use
%       a final image at step 4 and yield PIV results at step 3.5, or
%       alternatively, setting gap -> 2 would use a final image at step 5
%       and yield PIV results at step 4.
%   samp_len: see piv_step() help for details
%   samp_spc: see piv_step() help for details
%   intr_len: see piv_step() help for details
%   num_pass: see piv_step() help for details
%   valid_radius: see piv_step() help for details
%   valid_max: see piv_step() help for details
%   valid_eps: see piv_step() help for details
%   min_frac_data: see piv_step() help for details
%   min_frac_overlap: see piv_step() help for details
%   notes: String, notes to be included in output MAT-file as a global
%       attribute. default = ''
% %

update_path('piv', 'util');

% set defaults
narginchk(13, 14);
if nargin < 14; notes = ''; end

% check for sane arguments (pass-through arguments are checked in subroutines)
validateattributes(output_file, {'char'}, {'vector'});
validateattributes(input_file, {'char'}, {'vector'});
validateattributes(step_range, {'numeric'}, {'vector', 'numel', 2});
validateattributes(gap, {'numeric'}, {'scalar', 'positive', 'integer'});
[~, ~, output_file_ext] = fileparts(output_file);
assert(strcmp('.mat', output_file_ext), 'Output file must be .mat');

% fail if output exists already
assert(exist(output_file, 'file') == 0, ...
    'Output file exists, either make space or choose another filename');

% create output file 
% note: populated progressively to preserve results in the event of an error
output = matfile(output_file, 'Writable', true);

% read in all required input data from file
% note: reading a single time step takes about the same amount of time as
%   reading the whole dataset, so we do not read a step-at-a time
% note: it is possible this may cause memory problems, but ~500 image
%   series fits in ~5 GB, so no reason to solve this problem now
[x_img, y_img, step_img, img, mask_img] = read_input_data(input_file);

% compute indices for image pairs for given step range and gap
if isnan(step_range(1))
    min_idx = 1;
else
    min_idx = find(step_img >= step_range(1), 1, 'first');
end
if isnan(step_range(2))
    max_idx = length(step_img);
else
    max_idx = find(step_img <= step_range(2), 1, 'last');
end
ini_idx = 1:length(step_img);
fin_idx = ini_idx + gap;
in_range = ini_idx >= min_idx & ini_idx <= max_idx ...
         & fin_idx >= min_idx & fin_idx <= max_idx;
ini_idx = ini_idx(in_range);
fin_idx = fin_idx(in_range);

% compute step coordinate vector
step_piv = step_img(ini_idx);  % assign results to initial time
output.step = step_piv;  % save step coordinate vector
num_step = length(step_piv);

% define metadata
meta = struct();
meta.notes = notes;
meta.version = get_version();
meta.image_file_md5 = md5_hash(input_file);
meta.args.input_file = input_file;
meta.args.output_file = output_file;
meta.args.step_range = step_range;
meta.args.gap = gap;
meta.args.samp_len = samp_len;
meta.args.samp_spc = samp_spc;
meta.args.intr_len = intr_len;
meta.args.num_pass = num_pass;
meta.args.valid_radius = valid_radius;
meta.args.valid_max = valid_max;
meta.args.valid_eps = valid_eps;
meta.args.min_frac_data = min_frac_data;
meta.args.min_frac_overlap = min_frac_overlap;

meta.x.name = 'x';
meta.x.long_name = 'horizontal position, regular grid';
meta.x.notes = 'coordinate axis';
meta.x.dimensions = {};
meta.x.units = 'meters';

meta.y.name = 'y';
meta.y.long_name = 'vertical position, regular grid';
meta.y.notes = 'coordinate axis';
meta.y.dimensions = {};
meta.y.units = '';

meta.step.name = 'step';
meta.step.long_name = 'step number';
meta.step.notes = 'coordinate axis';
meta.step.dimensions = {};
meta.step.units = '1';

meta.u.name = 'u';
meta.u.long_name = 'displacement vector, x-component, regular grid';
meta.u.notes = '';
meta.u.dimensions = {'y', 'x', 'step'};
meta.u.units = 'meters/step';

meta.v.name = 'v';
meta.v.long_name = 'displacement vector, y-component, regular grid';
meta.v.notes = '';
meta.v.dimensions = {'y', 'x', 'step'};
meta.v.units = 'meters/step';

% note: build quality notes string dynamically from enumerated values
[qual_members, qual_names] = enumeration('Quality');
qual_values = Quality.to_uint8(qual_members);
qual_note = sprintf(['labels indicate how each observation was ', ...
    'computed, all values except %d (%s) are interpolated. Possible ', ...
    'values are: '], Quality.to_uint8(Quality.Valid), 'Valid');
for ii = 1:numel(qual_members)
    qual_note = [qual_note, sprintf('%d (%s), ', qual_values(ii), qual_names{ii})];  %#ok!
end
qual_note = sprintf('%s\b\b', qual_note);

meta.quality.name = 'quality';
meta.quality.long_name = 'quality labels';
meta.quality.notes = qual_note;
meta.quality.dimensions = {'x', 'y', 'step'};
meta.quality.units = 'categorical';

meta.mask.name = 'mask';
meta.mask.long_name = 'majority-sand mask';
meta.mask.notes = 'true where sample window contains majority-sand, false elsewhere';
meta.mask.dimensions = {'x', 'y', 'step'};
meta.quality.units = 'boolean';

output.meta = meta; % save to output file

% analyse all steps
for ii = 1:num_step
    
    fprintf('\n%s: begin step = %.1f\n', mfilename, step_piv(ii));
    fprintf('%s: ini_step = %d\n', mfilename, step_img(ini_idx(ii)));
    fprintf('%s: fin_step = %d\n', mfilename, step_img(fin_idx(ii)));
    
    % update image and mask pair
    img_ini = img(:, :, :, ini_idx(ii));
    mask_img_ini = mask_img(:, :, ini_idx(ii));
    
    img_fin = img(:, :, :, fin_idx(ii));
    mask_img_fin = mask_img(:, :, fin_idx(ii));
    
    % perform piv analysis
    piv_data = piv_step(...
        img_ini, img_fin, mask_img_ini, mask_img_fin, x_img, y_img, ...
        samp_len, samp_spc, intr_len, num_pass, valid_radius, ...
        valid_max, valid_eps, min_frac_data, min_frac_overlap); 
    
    if ii == 1
        % copy and save spatial coordinate vectors, created during PIV
        output.x = reshape(piv_data.x(1, :), size(piv_data.x, 2), 1);
        output.y = reshape(piv_data.y(:, 1), size(piv_data.y, 1), 1);
        % allocate other output variables
        % note: num_step is forced to be >= 2 b/c there is no way to have
        %   a singleton as the last dimension
        dims = [length(output.y), length(output.x), min(num_step, 2)];
        allocate(output, 'u', 'double', dims);
        allocate(output, 'v', 'double', dims);
        allocate(output, 'quality', 'uint8', dims);
        allocate(output, 'mask', 'logical', dims);
    end
    
    % save results for this step
    output.u(:, :, ii) = piv_data.u;
    output.v(:, :, ii) = piv_data.v;
    output.quality(:, :, ii) = Quality.to_uint8(piv_data.quality);
    output.mask(:, :, ii) = piv_data.mask;
    
    fprintf('%s: end step = %.1f\n', mfilename,  step_piv(ii));
end
fprintf('%s: end series\n', mfilename);


function [x, y, step, img, mask] = read_input_data(input_file)
% function [x, y, step, img, mask] = read_input_data(input_file)
% 
% Read in all required input data from file. Written as a function to scope
% memory use
% %
fprintf('%s: read input data from %s\n', mfilename, input_file);
input_data = load(input_file, 'x', 'y', 'step', 'img', 'mask');
x = double(input_data.x);
y = double(input_data.y);
step = double(input_data.step);
img = input_data.img;
mask = input_data.mask;
