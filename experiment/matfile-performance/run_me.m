% Setup a minimal test case that mimics the work we do in prep_series,
%   and then try a few different ways to solve it

NUM_STEPS = 508;
gcp();

% CASE_NAME = 'all-in-memory-parallel';
% + runs out of memory, just as expected

CASE_NAME = 'all-in-memory-serial';
% + memory goes to >90%
% + remarkably slow for a random number generator!




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
       
    % % write out to file at once
    % % note: skipping terrifically slow write step
    % fobj = matfile(output_file, 'Writable', true);
    % fobj.imgs = imgs;
    % fobj.masks = masks;
    % fobj.imgs_ext = imgs_ext;
    % fobj.masks_ext = masks_ext;
    
    end_case(CASE_NAME);
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
    rgb = randi(255, [nr, nc, 3], 'uint8');
end


function mask = make_mask(which)
    [nr, nc] = get_size(which);
    mask = logical(randi(1, [nr, nc]));
end


function eql = make_eql(which)
    [nr, nc] = get_size(which);
    eql = rand(nr, nc, 'single');
end