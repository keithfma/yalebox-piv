function [] = post_series(input_file, output_file, pro_bbox, ...
                          retro_bbox, pad_method, notes)
% function [] = post_series(input_file, output_file, pro_bbox, ...
%                           retro_bbox, pad_method, notes)
% 
% Run post-processing analyses on PIV data, and saves the results and
% metadata in a MAT-file file.
%
% Arguments:
%   input_file: string, file name of MAT-file containing PIV results as
%       produced by piv_series.m
% 
%   output_file: string, file name of MAT-file to create for post-processing
%       results data and metadata. Overwriting is not allowed to prevent
%       unhappy accidents
%
%   pro_bbox, retro_bbox: 4-element position vectors [xmin, ymin, width,
%       height] for the bounding boxes used to estimate pro- and retro-
%       plate displacements from PIV data.
%
%   pad_method: String, ...
%
%   notes: String, notes to be included in output MAT-file as a global
%       attribute. default = ''
% %

update_path('post', 'util');

% set defaults
if nargin < 6
    notes = '';
end

% check (immediate) input arguments
validateattributes(input_file, {'char'}, {'nonempty'}, ...
    mfilename, 'input_file');
assert(exist(input_file, 'file') == 2, ...
    sprintf('input file %s does not exist', input_file));
validateattributes(output_file, {'char'}, {'nonempty'}, ...
    mfilename, 'output_file');
assert(exist(output_file, 'file') ~= 2, ...
    sprintf('output file %s already exists', output_file));
validateattributes(notes, {'char'}, {}, mfilename, 'notes');

% create output file, fail if exists
assert(exist(output_file, 'file') == 0, ...
    'Output file exists, either make space or choose another filename');
output_data = matfile(output_file, 'Writable', true);

% define metadata
meta = struct();
meta.notes = notes;
meta.version = get_version();
meta.piv_file_md5 = md5_hash(input_file);
meta.args.pro_bbox = pro_bbox;
meta.args.retro_bbox = retro_bbox;
meta.args.pad_method = pad_method;

meta.x.name = 'x';
meta.x.long_name = 'horizontal position';
meta.x.notes = 'coordinate axis';
meta.x.dimensions = {};
meta.x.units = 'meters';

meta.y.name = 'y';
meta.y.long_name = 'vertical position';
meta.y.notes = 'coordinate axis';
meta.y.dimensions = {};
meta.y.units = 'meters';

meta.step.name = 'step';
meta.step.long_name = 'step_number';
meta.step.notes = 'coordinate axis';
meta.step.dimensions = {};
meta.step.units = '1';

meta.u_pro.name = 'u_pro';
meta.u_pro.long_name = 'median proside section displacement vector, x-component';
meta.u_pro.notes = 'sample from pixels inside pro_bbox';
meta.u_pro.dimensions = {'step'};
meta.u_pro.units = 'meters/step';

meta.v_pro.name = 'v_pro';
meta.v_pro.long_name = 'median proside section displacement vector, y-component';
meta.v_pro.notes = 'sample from pixels inside pro_bbox';
meta.v_pro.dimensions = {'step'};
meta.v_pro.units = 'meters/step';

meta.u_retro.name = 'u_retro';
meta.u_retro.long_name = 'median retroside section displacement vector, x-component';
meta.u_retro.notes = 'sample from pixels inside retro_bbox';
meta.u_retro.dimensions = {'step'};
meta.u_retro.units = 'meters/step';

meta.v_retro.name = 'v_retro';
meta.v_retro.long_name = 'median retroside section displacement vector, y-component';
meta.v_retro.notes = 'sample from pixels inside retro_bbox';
meta.v_retro.dimensions = {'step'};
meta.v_retro.units = 'meters/step';

meta.F11.name = 'F11';
meta.F11.long_name = 'dx/dX, deformation gradient tensor element (1,1)';
meta.F11.notes = '';
meta.F11.dimensions = {'y', 'x', 'step'};
meta.F11.units = '1';

meta.F12.name = 'F12';
meta.F12.long_name = 'dx/dY, deformation gradient tensor element (1,2)';
meta.F12.notes = '';
meta.F12.dimensions = {'y', 'x', 'step'};
meta.F12.units = '1';

meta.F21.name = 'F21';
meta.F21.long_name = 'dy/dX, deformation gradient tensor element (2,1)';
meta.F21.notes = '';
meta.F21.dimensions = {'y', 'x', 'step'};
meta.F21.units = '1';

meta.F22.name = 'F22';
meta.F22.long_name = 'dy/dY, deformation gradient tensor element (2,2)';
meta.F22.notes = '';
meta.F22.dimensions = {'y', 'x', 'step'};
meta.F22.units = '1';

meta.S1.name = 'S1';
meta.S1.long_name = 'maximum principle stretch';
meta.S1.notes = '';
meta.S1.dimensions = {'y', 'x', 'step'};
meta.S1.units = '1';

meta.S1x.name = 'S1x';
meta.S1x.long_name = 'maximum principle stretch unit direction vector, x-component';
meta.S1x.notes = '';
meta.S1x.dimensions = {'y', 'x', 'step'};
meta.S1x.units = '1';

meta.S1y.name = 'S1y';
meta.S1y.long_name = 'maximum principle stretch unit direction vector, y-component';
meta.S1y.notes = '';
meta.S1y.dimensions = {'y', 'x', 'step'};
meta.S1y.units = '1';

meta.S2.name = 'S2';
meta.S2.long_name = 'minimum principle stretch';
meta.S2.notes = '';
meta.S2.dimensions = {'y', 'x', 'step'};
meta.S2.units = '1';

meta.S2x.name = 'S2x';
meta.S2x.long_name = 'minimum principle stretch unit direction vector, x-component';
meta.S2x.notes = '';
meta.S2x.dimensions = {'y', 'x', 'step'};
meta.S2x.units = '1';

meta.S2y.name = 'S2y';
meta.S2y.long_name = 'minimum principle stretch unit direction vector, y-component';
meta.S2y.notes = '';
meta.S2y.dimensions = {'y', 'x', 'step'};
meta.S2y.units = '1';

meta.spin.name = 'spin';
meta.spin.long_name = 'local rotation angle, clockwise';
meta.spin.notes = '';
meta.spin.dimensions = {'y', 'x', 'step'};
meta.spin.units = 'radians';

meta.D1.name = 'D1';
meta.D1.long_name = 'equivalent maximum principle stretch rate';
meta.D1.notes = '';
meta.D1.dimensions = {'y', 'x', 'step'};
meta.D1.units = '1/step';

meta.D2.name = 'D2';
meta.D2.long_name = 'equivalent minimum principle stretch rate';
meta.D2.notes = '';
meta.D2.dimensions = {'y', 'x', 'step'};
meta.D2.units = '1/step';

meta.Dd.name = 'Dd';
meta.Dd.long_name = 'deviatoric strain rate';
meta.Dd.notes = 'see Brandon (1995), Ring and Brandon (1999)';
meta.Dd.dimensions = {'y', 'x', 'step'};
meta.Dd.units = '1/step';

meta.Dv.name = 'Dv';
meta.Dv.long_name = 'volume strain rate';
meta.Dv.notes = 'see Brandon (1995), Ring and Brandon (1999)';
meta.Dv.dimensions = {'y', 'x', 'step'};
meta.Dv.units = '1/step';

meta.Wk_star.name = 'Wk_star';
meta.Wk_star.long_name = 'kinematic vorticity number, relative to Dd';
meta.Wk_star.notes = 'see Brandon (1995), Ring and Brandon (1999)';
meta.Wk_star.dimensions = {'y', 'x', 'step'};
meta.Wk_star.units = '1';

meta.Ak_star.name = 'Ak_star';
meta.Ak_star.long_name = 'kinematic dilatancy number, relative to Dd';
meta.Ak_star.notes = 'see Brandon (1995), Ring and Brandon (1999)';
meta.Ak_star.dimensions = {'y', 'x', 'step'};
meta.Ak_star.units = '1';

output_data.meta = meta;

% open input file for reading, read in coordinate data
input_data = matfile(input_file, 'Writable', false);
xx   = double(input_data.x_grd); 
yy   = double(input_data.y_grd); 
step = double(input_data.step);
num_steps = length(step);
num_x = length(xx);
num_y = length(yy);

% write constant variables to output
output_data.x = xx;
output_data.y = yy;
output_data.step = step;

% allocate remaining output variables
vars_1d = {'u_pro', 'v_pro', 'u_retro', 'v_retro'};
for ii = 1:length(vars_1d)
    allocate(output_data, vars_1d{ii}, 'double', [num_steps, 1]);
end
vars_3d = {'F11', 'F12', 'F21', 'F22', 'S1', 'S1x', 'S1y', 'S2', ...
    'S2x', 'S2y', 'spin', 'D1', 'D2', 'Dd', 'Dv', 'Wk_star', 'Ak_star'};
for ii = 1:length(vars_3d)
    allocate(output_data, vars_3d{ii}, 'double', [num_y, num_x, num_steps]);
end

% run analyses for each timestep and save results
for ii = 1:num_steps
    
    uu = input_data.u_grd(:,:,ii);
    vv = input_data.v_grd(:,:,ii);
    roi = input_data.roi_grd(:,:,ii);
    
    uv_pro   = post_displ_rect(xx, yy, uu, vv, pro_bbox);
    output_data.u_pro(ii, 1) = uv_pro(1);
    output_data.v_pro(ii, 1) = uv_pro(2);
    
    uv_retro = post_displ_rect(xx, yy, uu, vv, retro_bbox);
    output_data.u_retro(ii, 1) = uv_retro(1);
    output_data.v_retro(ii, 1) = uv_retro(2);
    
    strain = post_strain(xx, yy, uu, vv, roi, pad_method);
    for jj = 1:length(vars_3d)
        output_data.(vars_3d{jj})(:, :, ii) = strain.(vars_3d{jj});
    end
    
    fprintf('%s: %d of %d\n', mfilename, ii, num_steps);
end



