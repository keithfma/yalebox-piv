%
% Run a set of experiments exploring a set of prep and/or piv parameter cases
%
% Arguments:
%   default_prep_param_file
%   default_piv_param_file
%   prep_case
%   piv_case
%   dest_dir
%   test_index
%   overwrite


% note: not bothering to include test cases for post-processing, as there are basically
%   no adjustable parameters at present.

update_path('util', 'prep', 'piv', 'post', 'allcomb');

PLOT_FORMATS = {'png', 'fig'};

% TEST: run as a script for now
default_prep_param_file = './experiment/choose-fill/prep-parameters.mat';
default_piv_param_file = './experiment/choose-fill/piv-parameters.mat';
prep_case = struct();
prep_case.fill_skin_min = 5;
prep_case.fill_skin_max = 10;
piv_case = struct();
index = 1;
dest_dir = './delete_me';
overwrite = false;
% END TEST

% apply cases
prep_param = update_params(default_prep_param_file, prep_case);
prep_param.exp_files.value = prep_param.exp_files.value(index:index+piv_param.gap.value);
piv_param = update_params(default_piv_param_file, piv_case);

% run prep
prep_name = get_case_name(prep_case);
image_file = get_file_name(dest_dir, prep_name, '', '-IMAGE.mat');
fpromptf('Prep: %s', prep_name);
if need_to_run(image_file, overwrite)
    prep(prep_param, image_file);
end

% run PIV
piv_name = get_case_name(piv_case);
piv_file = get_file_name(dest_dir, prep_name, piv_name, '-VELOCITY.mat');
fpromptf('PIV: %s', piv_name);
if need_to_run(piv_file, overwrite)
    piv(piv_param, image_file, piv_file);
end
piv_result = get_single_piv_result(load(piv_file));

% run strain
this_strain_result = post_strain(piv_result.x, piv_result.y, piv_result.u, piv_result.v);

% make plots
plot_prefix = sprintf('PIV: %s; PREP: %s;', piv_name, prep_name);
display_piv(piv_result, plot_prefix)
for kk = 1:length(PLOT_FORMATS)
    fn = get_file_name(dest_dir, prep_name, piv_name, sprintf('-VELOCITY.%s', PLOT_FORMATS{kk}));
    saveas(gcf, fn);
end

display_strain(this_strain_result, piv_result.mask, plot_prefix);
for kk = 1:length(PLOT_FORMATS)
    fn = get_file_name(dest_dir, prep_name, piv_name, sprintf('-STRAIN.%s', PLOT_FORMATS{kk}));
    saveas(gcf, fn);
end


function fpromptf(format, varargin)
% helper function for writing a big fat prompt

    fprintf('\n----------------------------------\n');
    fprintf([format, '\n'], varargin{:});
    fprintf('----------------------------------\n\n');

end


function name = get_file_name(folder, prep_name, piv_name, suffix)
% helper function to build a standardized filename
    
name = '';
    if ~isempty(prep_name)
        name = sprintf('PREP-%s', prep_name);
    end
    if ~isempty(piv_name)
        name = [sprintf('PIV-%s-', piv_name), name];
    end
    name = fullfile(folder, [name, suffix]);
    
end
    

function tf = need_to_run(file, overwrite)
% helper function to prepare for running a arbitrary stage, return True to run, or False to skip
     if isfile(file) && overwrite
            delete(file);
     end
     if isfile(file)
         fprintf('File exists, skipping: %s\n', file);
         tf = false;
     else
         tf = true;
     end
end
   

function result = get_single_piv_result(results)
% helper function to pull a single PIV result struct out of a series of them
    
    result = struct('x', results.x, 'y', results.y, 'u', results.u(:, :, 1), ...
                    'v', results.v(:, :, 1), 'mask', results.mask(:, :, 1), ...
                    'quality', results.quality(:, :, 1));

end


function params = update_params(base_param_file, case_params)
% helper function to update base parameters for each experiment case, taking care that the
%   parameters updated exist.

    params = load_param(base_param_file);
    
    case_param_names = fieldnames(case_params);
    for ii = 1:length(case_param_names)
        ensure_has_field(params, case_param_names{ii});
        ensure_has_field(params.(case_param_names{ii}), 'value');
        params.(case_param_names{ii}).value = case_params.(case_param_names{ii});
    end
    
end


function ensure_has_field(a_struct, a_field)
% helper function that fails if the input struct does not have the input field defined
    matched = cellfun(@(x) strcmp(x, a_field), fieldnames(a_struct), 'UniformOutput', true);
    try
        assert(sum(matched) == 1);
    catch
        keyboard
    end
end

    

function case_name = get_case_name(params)
% helper function to build a useful name out of a set of input parameters

    if isempty(fieldnames(params))
        case_name = 'default';
    else
        case_name = '';
        param_names = fieldnames(params);
        for ii = 1:length(param_names)
            case_name = extend_string(case_name, param_names{ii}, params.(param_names{ii}));
        end
    end
    
end


function str = extend_string(str, name, value)
% helper function for building up case names one peice at a time

    new = sprintf(['%s=', fmt(value)], name, value);
    if isempty(str)
        str = new;
    else
        str = [str, ',', new];
    end

end


function fmt_str = fmt(value)
% helper function for sprintf-style formatting when we do not know the data type

    if isnumeric(value) && isreal(value) && rem(value, 1) == 0
        fmt_str = '%d';
    elseif isfloat(value)
        fmt_str = '%.3f';
    elseif ischar(value)
        fmt_str = '%s';
    else
        error('What is this?');
    end

end