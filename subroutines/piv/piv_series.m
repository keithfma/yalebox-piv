function [] = piv_series(...
    output_file, input_file, step_range, gap, samp_len, samp_spc, ...
    intr_len, num_pass, valid_radius, valid_max, valid_eps, ...
    min_frac_data, min_frac_overlap, pad_method, notes)
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
% 
%   input_file: String, name of the input MAT file containing pre-processed
%       image data
%
%   step_range: 2-element vector, [initial, final] steps of input images series
%       to include in PIV analysis. Either element can be set to NaN to use the
%       full range.
%
%   gap : Scalar, integer, gap between analyzed images in steps, for
%       example, for an initial image at step 3, setting gap -> 1 would use
%       a final image at step 4 and yield PIV results at step 3.5, or
%       alternatively, setting gap -> 2 would use a final image at step 5
%       and yield PIV results at step 4.
%
%   samp_len: see piv_step() help for details
%
%   samp_spc: see piv_step() help for details
%
%   intr_len: see piv_step() help for details
%
%   num_pass: see piv_step() help for details
%
%   valid_radius: see piv_step() help for details
%
%   valid_max: see piv_step() help for details
%
%   valid_eps: see piv_step() help for details
%
%   min_frac_data: see piv_step() help for details
%
%   min_frac_overlap: see piv_step() help for details
%
%   pad_method: see piv_step() help for details
%
%   notes: String, notes to be included in output MAT-file as a global
%       attribute. default = ''
% %

update_path('piv', 'util');

% set defaults
narginchk(14, 15);
if nargin < 15; notes = ''; end

% check for sane arguments (pass-through arguments are checked in subroutines)
validateattributes(output_file, {'char'}, {'vector'});
validateattributes(input_file, {'char'}, {'vector'});
validateattributes(step_range, {'numeric'}, {'vector', 'numel', 2});
validateattributes(gap, {'numeric'}, {'scalar', 'positive', 'integer'});
[~, ~, output_file_ext] = fileparts(output_file);
assert(strcmp('.mat', output_file_ext), 'Output file must be .mat');

% open the input data file for reading
input_data = matfile(input_file, 'Writable', false);

% read constants (input coordinate vectors)
x_img = double(input_data.x);
y_img = double(input_data.y);
step_img = double(input_data.step);

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
step = step_img(ini_idx);
num_step = length(step);

% create output file, fail if exists
assert(exist(output_file, 'file') == 0, ...
    'Output file exists, either make space or choose another filename');
output_data = matfile(output_file, 'Writable', true);

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
meta.args.pad_method = pad_method;

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

output_data.meta = meta;

% analyse all steps
for ii = 1:num_step
    
    fprintf('\n%s: begin step = %.1f\n', mfilename, step(ii));
    fprintf('%s: ini_step = %d\n', mfilename, step_img(ini_idx(ii)));
    fprintf('%s: fin_step = %d\n', mfilename, step_img(fin_idx(ii)));
    
    % update image and mask pair
    img0 = input_data.img(:, :, :, ini_idx(ii));
    mask0 = input_data.mask(:, :, ini_idx(ii));
    
    img1 = input_data.img(:, :, :, fin_idx(ii));
    mask1 = input_data.mask(:, :, fin_idx(ii));
    
    % perform piv analysis
    piv_data = piv_step(img0, img1, mask0, mask1, x_img, y_img, samp_len, samp_spc, ...
        intr_len, num_pass, valid_radius, valid_max, valid_eps, ...
        min_frac_data, min_frac_overlap, pad_method); 
    
    % write results to output file
    if ii == 1
        % write spatial coordinate vectors,  created during PIV
        output_data.x = reshape(piv_data.x(1, :), size(piv_data.x, 2), 1);
        output_data.y = reshape(piv_data.y(:, 1), size(piv_data.y, 1), 1);
        output_data.step = step(:);        
        % allocate other output variables
        % note: add one to num_step to handle special case for single step test
        %   runs, as there is no way to have a singleton as the last dimension
        dimensions = [length(output_data.y), length(output_data.x), num_step+1];
        allocate(output_data, 'u', 'double', dimensions);
        allocate(output_data, 'v', 'double', dimensions);
        allocate(output_data, 'quality', 'uint8', dimensions);
    end
    output_data.u(:, :, ii) = piv_data.u;
    output_data.v(:, :, ii) = piv_data.v;
    output_data.quality(:, :, ii) = Quality.to_uint8(piv_data.quality);
    
    fprintf('%s: end step = %.1f\n', mfilename,  step(ii));
end

fprintf('%s: end series\n', mfilename);
fprintf('%s: output data saved to %s\n', mfilename, output_file);