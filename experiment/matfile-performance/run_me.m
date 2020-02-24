% Setup a minimal test case that mimics the work we do in prep_series,
%   and then try a few different ways to solve it

% add dependencies to path
addpath('../../dependencies/newmatic');

% NUM_STEPS = 508;
NUM_STEPS = 200;


% note: open the pool once, at the start of the script
gcp();

% note: random number generation is slow at this scale, but I think it is important
%   to randomize the array contents so that compression is not too easy and effective
rng(0, 'simdTwister');  % 2x faster than default

% CASE_NAME = 'all-in-memory-parallel';
% + runs out of memory, just as expected

% CASE_NAME = 'all-in-memory-serial';
% % + memory goes to ~100%

% CASE_NAME = 'array-serial';
% % + partial write basically causes the program to stall

% CASE_NAME = 'cell-array-serial';
% % + fails, docs say: "The MAT-file object does not support indexing into: ... Cells of cell arrays

% CASE_NAME = 'array-serial-rows';
% % + not better

CASE_NAME = 'array-serial-newmatic';
% + no memory problems
% + about 2x slower than in memory, but acceptable

if strcmp(CASE_NAME, 'all-in-memory-parallel') % --------------------------------------------------
    
    output_file = results_file(CASE_NAME);
    
    start_case(CASE_NAME);
    
    % init arrays
    [nr, nc] = get_size('orig');
    imgs = zeros(nr, nc, 3, NUM_STEPS, 'uint8');
    masks = false(nr, nc, NUM_STEPS);
    
    % populate arrays in parallel
    parfor ii = 1:NUM_STEPS
        imgs(:, :, :, ii) = make_rgb('orig');
        masks(:, :, ii) = make_mask('orig');
    end
    
    % trim the arrays to data extent
    [nr, nc] = get_size('trim');
    imgs = imgs(1:nr, 1:nc, :, :);
    masks = masks(1:nr, 1:nc, :);
    
    % init padded arrays
    [nr, nc] = get_size('pad');
    imgs_ext = zeros(nr, nc, NUM_STEPS, 'single');
    masks_ext = false(nr, nc, NUM_STEPS);
    
    % populate padded arrays in parallel
    parfor ii = 1:NUM_STEPS
        imgs_ext(:, :, ii) = make_eql('pad');
        masks_ext(:, :, ii) = make_mask('pad');
    end
       
    % write out to file at once
    fobj = matfile(output_file, 'Writable', true);
    fobj.imgs = imgs;
    fobj.masks = masks;
    fobj.imgs_ext = imgs_ext;
    fobj.masks_ext = masks_ext;
    
    end_case(CASE_NAME);

    
    
elseif strcmp(CASE_NAME, 'all-in-memory-serial') % ----------------------------------------------
    
    
    output_file = results_file(CASE_NAME);
    
    start_case(CASE_NAME);
    
    % init arrays
    [nr, nc] = get_size('orig');
    imgs = zeros(nr, nc, 3, NUM_STEPS, 'uint8');
    masks = false(nr, nc, NUM_STEPS);
    
    % populate arrays
    for ii = 1:NUM_STEPS
        imgs(:, :, :, ii) = make_rgb('orig');
        masks(:, :, ii) = make_mask('orig');
    end
    
    % trim the arrays to data extent
    [nr, nc] = get_size('trim');
    imgs = imgs(1:nr, 1:nc, :, :);
    masks = masks(1:nr, 1:nc, :);
    
    % init padded arrays
    [nr, nc] = get_size('pad');
    imgs_ext = zeros(nr, nc, NUM_STEPS, 'single');
    masks_ext = false(nr, nc, NUM_STEPS);
    
    % populate padded arrays
    for ii = 1:NUM_STEPS
        imgs_ext(:, :, ii) = make_eql('pad');
        masks_ext(:, :, ii) = make_mask('pad');
    end
       
    % write out to file at once
    % note: skipping terrifically slow write step
    fobj = matfile(output_file, 'Writable', true);
    fobj.imgs = imgs;
    fobj.masks = masks;
    fobj.imgs_ext = imgs_ext;
    fobj.masks_ext = masks_ext;
    
    end_case(CASE_NAME);
    
    
elseif strcmp(CASE_NAME, 'array-serial') % ----------------------------------------------
    
    
    temp_file = matfile(results_file('delete_me'), 'Writable', true);
    output_file = matfile(results_file(CASE_NAME), 'Writable', true);
    
    start_case(CASE_NAME);
    
    % init temp file
    [nr, nc] = get_size('orig');
    allocate(temp_file, 'imgs', 'uint8', [nr, nc, 3, NUM_STEPS]);
    allocate(temp_file, 'masks', 'logical', [nr, nc, NUM_STEPS]);
    
    % populate temp file
    for ii = 1:NUM_STEPS
        temp_file.imgs(:, :, :, ii) = make_rgb('orig');
        temp_file.masks(:, :, ii) = make_mask('orig');
    end
    
    % init results file
    [nr, nc] = get_size('trim');
    allocate(output_file, 'imgs', 'uint8', [nr, nc, 3, NUM_STEPS]);
    allocate(output_file, 'masks', 'logical', [nr, nc, NUM_STEPS]);
    
    [nr, nc] = get_size('pad');
    allocate(output_file, 'imgs_ext', 'single', [nr, nc, NUM_STEPS]);
    allocate(output_file, 'masks_ext', 'logical', [nr, nc, NUM_STEPS]);
    
    % populate results file
    % note: not really fair, since we should be reading from file at each step here,
    %   but that doen't seem like it is worth the complexity of testing...
    for ii = 1:NUM_STEPS
        output_file.imgs(:, :, :, ii) = make_rgb('trim');
        output_file.masks(:, :, ii) = make_mask('trim');
        output_file.imgs_ext(:, :, ii) = make_eql('pad');
        output_file.masks_ext(:, :, ii) = make_mask('pad');
    end
    
    end_case(CASE_NAME);

    
elseif strcmp(CASE_NAME, 'cell-array-serial') % ----------------------------------------------
    
    
    temp_file = matfile(results_file('delete_me'), 'Writable', true);
    output_file = matfile(results_file(CASE_NAME), 'Writable', true);
    
    start_case(CASE_NAME);
    
    % init temp file
    temp_file.imgs = cell(NUM_STEPS, 1);
    temp_file.masks = cell(NUM_STEPS, 1);
    
    % populate temp file
    for ii = 1:NUM_STEPS
        temp_file.imgs(ii, 1) = make_rgb('orig');
        temp_file.masks(ii, 1) = make_mask('orig');
    end
    
    % init results file
    [nr, nc] = get_size('trim');
    output_file.imgs = cell(NUM_STEPS, 1);
    output_file.masks = cell(NUM_STEPS, 1);
    output_file.imgs_ext = cell(NUM_STEPS, 1);
    output_file.masks_ext = cell(NUM_STEPS, 1);
    
    % populate results file
    % note: not really fair, since we should be reading from file at each step here,
    %   but that doen't seem like it is worth the complexity of testing...
    for ii = 1:NUM_STEPS
        output_file.imgs(ii, 1) = make_rgb('trim');
        output_file.masks(ii, 1) = make_mask('trim');
        output_file.imgs_ext(ii, 1) = make_eql('pad');
        output_file.masks_ext(ii, 1) = make_mask('pad');
    end
       
    end_case(CASE_NAME);

    
elseif strcmp(CASE_NAME, 'array-serial-rows') % ----------------------------------------------
    
    
    temp_file = matfile(results_file('delete_me'), 'Writable', true);
    output_file = matfile(results_file(CASE_NAME), 'Writable', true);
    
    start_case(CASE_NAME);
    
    % init temp file
    [nr, nc] = get_size('orig');
    allocate(temp_file, 'imgs', 'uint8', [NUM_STEPS, nr, nc, 3]);
    allocate(temp_file, 'masks', 'logical', [NUM_STEPS, nr, nc]);
    
    % populate temp file
    [nr, nc] = get_size('orig');
    for ii = 1:NUM_STEPS
        temp_file.imgs(ii, :, :, :) = reshape(make_rgb('orig'), [1, nr, nc, 3]);
        temp_file.masks(ii, :, :) = reshape(make_mask('orig'), [1, nr, nc]);
    end
    
    % init results file
    [nr, nc] = get_size('trim');
    allocate(output_file, 'imgs', 'uint8', [NUM_STEPS, nr, nc, 3]);
    allocate(output_file, 'masks', 'logical', [NUM_STEPS, nr, nc]);
    
    [nr, nc] = get_size('pad');
    allocate(output_file, 'imgs_ext', 'single', [NUM_STEPS, nr, nc]);
    allocate(output_file, 'masks_ext', 'logical', [NUM_STEPS, nr, nc]);
    
    % populate results file
    % note: not really fair, since we should be reading from file at each step here,
    %   but that doen't seem like it is worth the complexity of testing...
    for ii = 1:NUM_STEPS
        output_file.imgs(ii, :, :, :) = make_rgb('trim');
        output_file.masks(ii, :, :) = make_mask('trim');
        output_file.imgs_ext(ii, :, :) = make_eql('pad');
        output_file.masks_ext(ii, :, :) = make_mask('pad');
    end
    
    end_case(CASE_NAME);

        
elseif strcmp(CASE_NAME, 'array-serial-newmatic') % ----------------------------------------------
    
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
        output_file.imgs(:, :, :, ii) = make_rgb('trim');
        output_file.masks(:, :, ii) = make_mask('trim');
        output_file.imgs_ext(:, :, ii) = make_eql('pad');
        output_file.masks_ext(:, :, ii) = make_mask('pad');
    end
    
    end_case(CASE_NAME);
   
else
    error('Unknown case name');
    
end


% --------------------------------------------------------------------------------------------------

   

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