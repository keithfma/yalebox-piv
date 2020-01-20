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

update_path('util', 'prep', 'piv', 'allcomb');

% TEST: cases as a struct
base_prep_param_file = './experiment/choose-fill/prep-parameters.mat';

prep_test_params = struct();
prep_test_params.fill_skin_min = {3, 5};
prep_test_params.fill_skin_max = {10, 15, 20};


dest_dir = './delete_me';

prep_cases = get_cases(prep_test_params);

prep_results = cell(size(prep_cases));
for ii = 1:length(prep_cases)
   
    this_case = prep_cases(ii);
    this_name = get_case_name(this_case);
    this_image_file = fullfile(dest_dir, sprintf('IMAGE-%s.mat', this_name));
    this_param = update_case_params(base_prep_param_file, this_case);
    
    fprintf('----------------------------------\n');
    fprintf('Prep case %d/%d: %s\n', ii, length(prep_cases), this_name);
    fprintf('----------------------------------\n\n');
    prep(this_param, this_image_file);
    prep_results{ii} = this_image_file;
end


% prep_cases_names = cell(size(prep_cases));
% for ii = 1:num_cases
%     prep_cases_names{ii} = get_case_name(prep_cases(ii));
% end





% So, we are getting closer here, but, how to deal with PAIRS of PIV and PREP params?
% + Seems like we need a param file for each case in PREP and PIV, 


% TODO: build a cell array of prep input files

% TODO: build a cell array of piv input files in the same way (generalize)

% TODO: run all prep cases

% TODO: for each prep case, run al PIV cases


function cases = get_cases(params)
% helper function to unpack all combinations of the test parameter definitions into
%   an array of structs with a single value for each parameter
    
    if isempty(params)
        % handle case where no parameters are to be tested for this step (prep or piv)
        cases = struct();
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

    

function case_name = get_case_name(params, prefix)
% helper function to build a useful name out of a set of input parameters

    case_name = '';
    param_names = fieldnames(params);
    for ii = 1:length(param_names)
        case_name = extend_string(case_name, param_names{ii}, params.(param_names{ii}));
    end
    if nargin == 2
        case_name = [prefix, case_name];
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