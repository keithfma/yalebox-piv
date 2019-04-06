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
%   samp_len: see piv() help for details
%
%   samp_spc: see piv() help for details
%
%   intr_len: see piv() help for details
%
%   num_pass: see piv() help for details
%
%   valid_radius: see piv() help for details
%
%   valid_max: see piv() help for details
%
%   valid_eps: see piv() help for details
%
%   min_frac_data: see piv() help for details
%
%   min_frac_overlap: see piv() help for details
%
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

% open the input data file for reading
input_data = matfile(input_file, 'Writable', false);

% read constants (input coordinate vectors)
x_img = double(input_data.x);
y_img = double(input_data.y);
step_img = double(input_data.step);

% compute indices for image pairs
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

% compute size of gridded output dimensions
[r_grd, c_grd] = piv_sample_grid(samp_len, samp_spc, length(y_img), length(x_img));
step = step_img(ini_idx);
num_x_grd = size(c_grd, 2);
num_y_grd = size(r_grd, 1);
num_step = length(step);
max_num_pts = num_x_grd*num_y_grd;

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

meta.x_grd.name = 'x_grd';
meta.x_grd.long_name = 'horizontal position, regular grid';
meta.x_grd.notes = 'coordinate axis';
meta.x_grd.dimensions = {};
meta.x_grd.units = 'meters';

meta.y_grd.name = 'y_grd';
meta.y_grd.long_name = 'vertical position, regular grid';
meta.y_grd.notes = 'coordinate axis';
meta.y_grd.dimensions = {};
meta.y_grd.units = '';

meta.step.name = 'step';
meta.step.long_name = 'step number';
meta.step.notes = 'coordinate axis';
meta.step.dimensions = {};
meta.step.units = '1';

meta.u_grd.name = 'u_grd';
meta.u_grd.long_name = 'displacement vector, x-component, interpolated to regular grid';
meta.u_grd.notes = '';
meta.u_grd.dimensions = {'y_grd', 'x_grd', 'step'};
meta.u_grd.units = 'meters/step';

meta.v_grd.name = 'v_grd';
meta.v_grd.long_name = 'displacement vector, y-component, interpolated to regular grid';
meta.v_grd.notes = '';
meta.v_grd.dimensions = {'y_grd', 'x_grd', 'step'};
meta.v_grd.units = 'meters/step';

meta.roi_grd.name = 'roi_grd';
meta.roi_grd.long_name = 'displacement vector mask, regular grid';
meta.roi_grd.notes = '';
meta.roi_grd.dimensions = {'y_grd', 'x_grd', 'step'};
meta.roi_grd.units = 'boolean';

meta.x_pts.name = 'x_pts';
meta.x_pts.long_name = 'horizontal position for raw measurements at scattered points';
meta.x_pts.notes = 'variable-length data stored as a fixed-size array, ignore NaNs';
meta.x_pts.dimensions = {'variable', 'step'};
meta.x_pts.units = 'meters';

meta.y_pts.name = 'y_pts';
meta.y_pts.long_name = 'vertical position for raw measurements at scattered points';
meta.y_pts.notes = 'variable-length data stored as a fixed-size array, ignore NaNs';
meta.y_pts.dimensions = {'variable', 'step'};
meta.y_pts.units = 'meters';

meta.u_pts.name = 'u_pts';
meta.u_pts.long_name = 'displacement vector, x-component, raw measurements at scattered points';
meta.u_pts.notes = 'variable-length data stored as a fixed-size array, ignore NaNs';
meta.u_pts.dimensions = {'variable', 'step'};
meta.u_pts.units = 'meters/step';

meta.v_pts.name = 'v_pts';
meta.v_pts.long_name = 'displacement vector, y-component, raw measurements at scattered points';
meta.v_pts.notes = 'variable-length data stored as a fixed-size array, ignore NaNs';
meta.v_pts.dimensions = {'variable', 'step'};
meta.v_pts.units = 'meters/step';
output_data.meta = meta;

% allocate output variables
% note: matfile does not allow indexing into cell arrays, so we have to
%   store the variable-length output in a fixed dimension grid
allocate(output_data, 'u_grd', 'double', [num_y_grd, num_x_grd, num_step]);
allocate(output_data, 'v_grd', 'double', [num_y_grd, num_x_grd, num_step]);
allocate(output_data, 'roi_grd', 'double', [num_y_grd, num_x_grd, num_step]);
allocate(output_data, 'x_pts', 'double', [max_num_pts, num_step]);
allocate(output_data, 'y_pts', 'double', [max_num_pts, num_step]);
allocate(output_data, 'u_pts', 'double', [max_num_pts, num_step]);
allocate(output_data, 'v_pts', 'double', [max_num_pts, num_step]);

% analyse all steps
for ii = 1:num_step
    
    fprintf('\n%s: begin step = %.1f\n', mfilename, step(ii));
    fprintf('%s: ini_step = %.1f\n', mfilename, step_img(ini_idx(ii)));
    fprintf('%s: fin_step = %.1f\n', mfilename, step_img(fin_idx(ii)));
    
    % update image and roi pair
    img0 = input_data.img(:, :, :, ini_idx(ii));
    roi0 = input_data.mask(:, :, ini_idx(ii));
    
    img1 = input_data.img(:, :, :, fin_idx(ii));
    roi1 = input_data.mask(:, :, fin_idx(ii));
    
    % perform piv analysis
    piv_data = piv_step(img0, img1, roi0, roi1, x_img, y_img, samp_len, samp_spc, ...
        intr_len, num_pass, valid_radius, valid_max, valid_eps, ...
        min_frac_data, min_frac_overlap); 
    
    % write results to output file
    % note: spatial coordinate vectors are created during PIV
    if ii == 1      
        output_data.x_grd = reshape(piv_data.x_grd(1, :), size(piv_data.x_grd, 2), 1);
        output_data.y_grd = reshape(piv_data.y_grd(:, 1), size(piv_data.y_grd, 1), 1);
        output_data.step = step(:);
    end
    output_data.u_grd(:, :, ii) = piv_data.u_grd;
    output_data.v_grd(:, :, ii) = piv_data.v_grd;
    
    fprintf('%s: end step = %.1f\n', mfilename,  step(ii));
end

fprintf('%s: end series\n', mfilename);
fprintf('%s: output data saved to %s\n', mfilename, output_file);