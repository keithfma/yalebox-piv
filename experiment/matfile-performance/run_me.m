% Setup a minimal test case that mimics the work we do in prep_series,
%   and then try a few different ways to solve it

% add dependencies to path
addpath('../../dependencies/newmatic');

% NUM_STEPS = 508;
NUM_STEPS = 30;
WAIT_TIME = 1;  % seconds, simulate a slower process for data generation

% note: open the pool once, at the start of the script
gcp();

% note: random number generation is slow at this scale, but I think it is important
%   to randomize the array contents so that compression is not too easy and effective
rng(0, 'simdTwister');  % 2x faster than default

% CASE_NAME = 'newmatic-serial';
% % + adding processing time makes this substantially slower

% CASE_NAME = 'afterEach';
% % + afterEach callbacks run serially, meaning we *can* use them for IO
% % + effectively hides processing time, this is the new winner

CASE_NAME = 'newmatic-parallel';
% + working...


if strcmp(CASE_NAME, 'newmatic-serial') % ----------------------------------------------
    
    start_case(CASE_NAME);
    
    % init temp_file
    [nr_orig, nc_orig] = get_size('orig');
    temp_file = newmatic(...
        results_file('delete_me'), ...
        newmatic_variable('imgs', 'uint8', [nr_orig, nc_orig, 3, NUM_STEPS], [nr_orig, nc_orig, 3, 1]), ...
        newmatic_variable('masks', 'logical', [nr_orig, nc_orig, NUM_STEPS], [nr_orig, nc_orig, 1]) ...
        );
    
    % populate temp file
    for ii = 1:NUM_STEPS
        pause(WAIT_TIME);  % pretend this takes a while
        temp_file.imgs(:, :, :, ii) = make_rgb('orig');
        temp_file.masks(:, :, ii) = make_mask('orig');
    end
    
    % init output file
    [nr_trim, nc_trim] = get_size('trim');
    [nr_pad, nc_pad] = get_size('pad');
    output_file = newmatic(...
        results_file(CASE_NAME), ...
        newmatic_variable('imgs', 'uint8', [nr_trim, nc_trim, 3, NUM_STEPS], [nr_trim, nc_trim, 3, 1]), ...
        newmatic_variable('masks', 'logical', [nr_trim, nc_trim, NUM_STEPS], [nr_trim, nc_trim, 1]), ...
        newmatic_variable('imgs_ext', 'single', [nr_pad, nc_pad, NUM_STEPS], [nr_pad, nc_pad, 1]), ...
        newmatic_variable('masks_ext', 'logical', [nr_pad, nc_pad, NUM_STEPS], [nr_pad, nc_pad, 1]) ...
        );
    
    % populate results file
    % note: not really fair, since we should be reading from file at each step here,
    %   but that doen't seem like it is worth the complexity of testing...
    for ii = 1:NUM_STEPS
        pause(WAIT_TIME);  % pretend this takes a while
        output_file.imgs(:, :, :, ii) = make_rgb('trim');
        output_file.masks(:, :, ii) = make_mask('trim');
        output_file.imgs_ext(:, :, ii) = make_eql('pad');
        output_file.masks_ext(:, :, ii) = make_mask('pad');
    end
    
    end_case(CASE_NAME);
    
elseif strcmp(CASE_NAME, 'afterEach')  % --------------------------------------------------------
    
    queue = parallel.pool.DataQueue;
    afterEach(queue, @sleep_and_print);
    
    parfor ii = 1:20
        send(queue, ii)
    end
    disp('loop complete');
    
    
elseif strcmp(CASE_NAME, 'newmatic-parallel') % ----------------------------------------------
    
    
    start_case(CASE_NAME);
    
    % init temp_file
    [nr_orig, nc_orig] = get_size('orig');
    temp_file = newmatic(...
        results_file('delete_me'), ...
        newmatic_variable('imgs', 'uint8', [nr_orig, nc_orig, 3, NUM_STEPS], [nr_orig, nc_orig, 3, 1]), ...
        newmatic_variable('masks', 'logical', [nr_orig, nc_orig, NUM_STEPS], [nr_orig, nc_orig, 1]) ...
    );
    
    % populate temp file
    temp_queue = parallel.pool.DataQueue;
    afterEach(temp_queue, @write_to_temp);
    
    parfor ii = 1:NUM_STEPS
        pause(WAIT_TIME);  % pretend this takes a while
        send(temp_queue, {temp_file, ii, make_rgb('orig'), make_mask('orig')});
    end
    
    % init output file
    [nr_trim, nc_trim] = get_size('trim');
    [nr_pad, nc_pad] = get_size('pad');
    output_file = newmatic(...
        results_file(CASE_NAME), ...
        newmatic_variable('imgs', 'uint8', [nr_trim, nc_trim, 3, NUM_STEPS], [nr_trim, nc_trim, 3, 1]), ...
        newmatic_variable('masks', 'logical', [nr_trim, nc_trim, NUM_STEPS], [nr_trim, nc_trim, 1]), ...
        newmatic_variable('imgs_ext', 'single', [nr_pad, nc_pad, NUM_STEPS], [nr_pad, nc_pad, 1]), ...
        newmatic_variable('masks_ext', 'logical', [nr_pad, nc_pad, NUM_STEPS], [nr_pad, nc_pad, 1]) ...
    );
    
    % populate results file
    output_queue = parallel.pool.DataQueue;
    afterEach(output_queue, @write_to_output);
    
    parfor ii = 1:NUM_STEPS
        pause(WAIT_TIME);  % pretend this takes a while
        send(output_queue, {...
            output_file, ii, make_rgb('trim'), make_mask('trim'), make_eql('pad'), make_mask('pad')});
    end
    
    end_case(CASE_NAME);
    
else
    
    error('Unknown case name');
    
end


% --------------------------------------------------------------------------------------------------


function write_to_temp(data)
    % data is {temp_file, ii, img, mask}
    temp_file = data{1};
    idx = data{2};
    temp_file.imgs(:, :, :, idx) = data{3};
    temp_file.masks(:, :, idx) = data{4};
    fprintf('temp %d\n', idx);
end


function write_to_output(data)
    % data is {output_file, ii, img, mask, img_ext, mask_ext}
    output_file = data{1};
    idx = data{2};
    output_file.imgs(:, :, :, idx) = data{3};
    output_file.masks(:, :, idx) = data{4};
    output_file.imgs_ext(:, :, idx) = data{5};
    output_file.masks_ext(:, :, idx) = data{6};
    fprintf('output %d\n', idx);
end


function sleep_and_print(val)
    pause(1);
    disp(val);
end



function start_case(name)
    fprintf('%s ------ start\n', name);
    tic
end


    function end_case(name)
        toc
        fprintf('%s ------ end\n', name);
    end


    function fn = results_file(name)
        % return path for results file, wiping out existing first
        fn = fullfile('.', 'results', sprintf('case-%s.mat', name));
        if isfile(fn)
            delete(fn)
        end
    end


    function [nr, nc] = get_size(which)
        switch which
            case 'orig'
                nr = 762;
                nc = 4277;
            case 'trim'
                nr = 500;
                nc = 4277;
            case 'pad'
                nr = 700;
                nc = 4477;
        end
    end


    function rgb = make_rgb(which)
        [nr, nc] = get_size(which);
        % rgb = randi(255, [nr, nc, 3], 'uint8');
        rgb = ones([nr, nc, 3], 'uint8');
    end


    function mask = make_mask(which)
        [nr, nc] = get_size(which);
        % mask = logical(randi(1, [nr, nc]));
        mask = true(nr, nc);
    end


    function eql = make_eql(which)
        [nr, nc] = get_size(which);
        % eql = rand(nr, nc, 'single');
        eql = ones(nr, nc, 'single');
    end