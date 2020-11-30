% Helper script for running PIV and subsequent processing. Purpose is to facilitate
%   twiddling the various parameters knobs and saving results with informative names.
% % 


+% note: assumes you ran update_path('main') in the root of yalebox_piv, so that the update_path
+%   function is available to run here
update_path('piv', 'util');

LOG_FILE = 'piv.log';
PIV_INPUT_PARAM_FILE_PREFIX = 'piv-input-param';


% define PIV input parameters 
p = piv_template();
p.notes.value = "Run PIV for Mark Brandon's prototype padding routine";
p.step_range.value = [100, 101];
p.samp_len.value = [100, 100];
p.intr_len.value = [200, 120];
p.num_pass.value = [1, 1];
p.samp_spc.value = 50;
p.valid_radius.value = 200;

% assign a unique case name
case_name = uuid(8);
write_params_to_log(LOG_FILE, case_name, p); 