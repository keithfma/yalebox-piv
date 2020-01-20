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

for ii = 1:length(prep_cases)
   
    % run prep case
    this_prep_case = prep_cases(ii);
    this_prep_name = get_case_name(this_prep_case);
    this_image_file = fullfile(dest_dir, sprintf('IMAGE-%s.mat', this_prep_name));
    this_prep_param = update_case_params(base_prep_param_file, this_prep_case);
    this_prep_param.exp_files.value = this_prep_param.exp_files.value(test_index:test_index+1);
    
    fprintf('----------------------------------\n');
    fprintf('Prep case %d/%d: %s\n', ii, length(prep_cases), this_prep_name);
    fprintf('----------------------------------\n\n');
    
    if isfile(this_image_file) && overwrite_prep
        delete(this_image_file)
    end
    if isfile(this_image_file)
        fprintf('Image file exists, skipping: %s\n', this_image_file);
    else
        prep(this_prep_param, this_image_file);
    end
    
    % TODO: plot any prep results you wish to plot
    
        
    for jj = 1:length(piv_cases)
        
        % run PIV case
        this_piv_case = piv_cases{jj};
        this_piv_name = get_case_name(this_piv_case);
        this_piv_file = fullfile(dest_dir, sprintf('PIV-%s-PREP-%s.mat', this_piv_name, this_prep_name));
        this_piv_param = update_case_params(base_piv_param_file, this_piv_case);
        
        fprintf('----------------------------------\n');
        fprintf('PIV case %d/%d: %s\n', jj, length(piv_cases), this_piv_name);
        fprintf('----------------------------------\n\n');
        
        if isfile(this_piv_file) && overwrite_piv
            delete(this_piv_file);
        end
        if isfile(this_piv_file)
            fprintf('PIV file exists, skipping: %s\n', this_piv_file);
        else
           piv(this_piv_param, this_image_file, this_piv_file);
        end
        
        % plot and save PIV results
        this_piv_result = get_single_piv_result(load(this_piv_file));
        this_plot_prefix = sprintf('PIV: %s; PREP: %s;', this_piv_name, this_prep_name);
        this_piv_plot_file = fullfile(...
            dest_dir, sprintf('PIV-%s-PREP-%s-DISPLACEMENT.png', this_piv_name, this_prep_name));
        display_piv(this_piv_result, this_plot_prefix)
        saveas(gcf, this_piv_plot_file, 'png');
        
        % TODO: loop over strain cases
        
        % TODO: use strain results saved above
        % plot and save strain results
        this_strain_plot_file = fullfile(...
            dest_dir, sprintf('PIV-%s-PREP-%s-STRAIN.png', this_piv_name, this_prep_name));
        this_strain_result = post_strain(...
            this_piv_result.x(1,:), ...
            this_piv_result.y(:,1), ...
            this_piv_result.u, ...
            this_piv_result.v);
        display_strain(this_strain_result, this_piv_result.mask, this_plot_prefix);
        saveas(gcf, this_strain_plot_file, 'png');
        
        keyboard
    end
    
end


function result = get_single_piv_result(results)
% helper function to pull a single PIV result struct out of a series of them
    
    [x, y] = meshgrid(results.x, results.y);
    result = struct('x', x, 'y', y, 'u', results.u(:, :, 1), ...
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