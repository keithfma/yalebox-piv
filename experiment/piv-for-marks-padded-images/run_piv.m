% Helper script for running PIV and subsequent processing. Purpose is to facilitate
%   twiddling the various parameters knobs and saving results with informative names.
% % 


% note: assumes you ran update_path('main') in the root of yalebox_piv, so that the update_path
%   function is available to run here
update_path('piv', 'util');

CASE_NAME = 'better-params';

CASE_FOLDER = sprintf('run_%s', CASE_NAME);
PIV_PARAM_FILE = sprintf('%s%spiv-param_%s.mat', CASE_FOLDER, filesep, CASE_NAME);
PIV_PARAM_LOG_FILE = sprintf('%s%spiv-param_%s.json', CASE_FOLDER, filesep, CASE_NAME);
PIV_DATA_FILE = sprintf('%s%spiv-data_%s.mat', CASE_FOLDER, filesep, CASE_NAME);
IMAGE_DATA_FILE = 'images.mat';


% define PIV input parameters 
p = piv_template();
p.notes.value = "Run PIV for Mark Brandon's prototype padding routine";
p.step_range.value = [100, 101];
p.samp_len.value = [100, 60, 30];
p.intr_len.value = [200, 120, 60];
p.num_pass.value = [1, 1, 1];
p.samp_spc.value = 15;
p.valid_radius.value = 45;

% clear and re-create output folder
mkdir(CASE_FOLDER);
save_param(p, PIV_PARAM_FILE);
param_to_json(PIV_PARAM_LOG_FILE, p, true);

% run!
piv(p, IMAGE_DATA_FILE, PIV_DATA_FILE);


