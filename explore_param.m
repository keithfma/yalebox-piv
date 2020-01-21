% function cases(base_prep_param_file, base_piv_param_file, prep_cases_map, piv_cases_map)
%
%
% Run a set of experiments exploring a set of prep and/or piv parameter cases
%
% Arguments:
%   base_prep_param_file
%   base_piv_param_file
%   prep_test_params
%   piv_test_params
%   dest_dir
%   test_index
%   overwrite_prep
%   overwrite_piv

% TODO: consider running from a struct, which I could save and load
% TODO: respect gap parameter
% TODO: remove a lot of this_ prefixes
% TODO: save all as fig and png (for exploration and sharing)
% TODO: drop allcomb, just run a single case and defer looping to outside

% note: not bothering to include test cases for post-processing, as there are basically
%   no adjustable parameters at present.

update_path('util', 'prep', 'piv', 'post', 'allcomb');

% TEST: run as a script for now
base_prep_param_file = './experiment/choose-fill/prep-parameters.mat';
base_piv_param_file = './experiment/choose-fill/piv-parameters.mat';

prep_test_params = struct();
prep_test_params.fill_skin_min = {3, 5};
prep_test_params.fill_skin_max = {10, 15, 20};

piv_test_params = struct();

test_index = 1;
dest_dir = './delete_me';
overwrite_prep = false;
overwrite_piv = false;
% END TEST

prep_cases = get_cases(prep_test_params);
piv_cases = get_cases(piv_test_params);

PLOT_FORMATS = {'png', 'fig'};

for ii = 1:length(prep_cases)
    
    % run prep
    prep_case = prep_cases(ii);
    prep_name = get_case_name(prep_case);
    image_file = get_file_name(dest_dir, prep_name, '', '-IMAGE.mat');
    
    fprintf('----------------------------------\n');
    fprintf('Prep case %d/%d: %s\n', ii, length(prep_cases), prep_name);
    fprintf('----------------------------------\n\n');
    
    if need_to_run(image_file, overwrite_prep)
        prep_param = update_case_params(base_prep_param_file, prep_case);
        prep_param.exp_files.value = prep_param.exp_files.value(test_index:test_index+1);
        prep(prep_param, image_file);
    end
        
    for jj = 1:length(piv_cases)
        
        % run PIV
        piv_case = piv_cases{jj};
        piv_name = get_case_name(piv_case);
        piv_file = get_file_name(dest_dir, prep_name, piv_name, '-VELOCITY.mat');
       
        fprintf('----------------------------------\n');
        fprintf('PIV case %d/%d: %s\n', jj, length(piv_cases), piv_name);
        fprintf('----------------------------------\n\n');
        
        if need_to_run(piv_file, overwrite_piv)
            piv_param = update_case_params(base_piv_param_file, piv_case);
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
           
    end
    
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


function cases = get_cases(params)
% helper function to unpack all combinations of the test parameter definitions into
%   an array of structs with a single value for each parameter
    
    if isempty(fieldnames(params))
        % handle case where no parameters are to be tested for this step (prep or piv)
        cases = cell(1, 1);
        cases{1} = struct();
        return
    end

    names = fieldnames(params);
    num_names = length(names);

    values = struct2cell(params);

    value_combos = allcomb(values{:});
    num_combos = length(value_combos);

    cases = cell(num_combos, 1);
    for ii = 1:num_combos
        this = struct();
        for jj = 1:num_names
            this.(names{jj}) = value_combos{ii, jj};
        end
        cases{ii} = this;
    end

    cases = cell2mat(cases);

end


function params = update_case_params(base_param_file, case_params)
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