function write_params_to_log(log_file, case_name, parameters)
% function write_params_to_log(log_file, case_name, parameters)
%
% Record case name and input parameters to log file, formatted as concatenated JSON objects.
% All parameters are included, but the .help text is dropped for brevity.
% 
% Arguments:
%   log_file: file name to append log entry to
%   case_name: a string name to attach to this log entry, typically, this would also
%       get attached to other outputs from whatever you are running (datafiles, figures, etc).
%   parameters: a yalebox-piv parameters struct (prep, piv, etc all should work)
% % 

% TODO: This does a lot more than we need it too! What we need is simpler, a function that
%   serializes and pretty-formats the input parameters, really, that is it. The case should 
%   probably be user defined, since hashes or uuids don't add much. 


% flatten (i.e., drop .help and promote .value)
flat_parameters = struct();

field_names = fieldnames(parameters);
for kk = 1:length(field_names)
    field_name = field_names{kk};
    flat_parameters.(field_name) = parameters.(field_name).value;
end

% nest parameters and case name side-by-side
out = struct( ...
    'case', case_name, ...
    'parameters', flat_parameters ...
);

% use built-in tools to serialize as JSON
raw_json = jsonencode(out, 'ConvertInfAndNaN', true);


% use external tools to format (specifically jq)
% note: annoyingly, we have to make a temporary file to pass JSON to `jq`, this is because MATLAB
%   passes everything directly to the shell, which then blows up if characters are not correctly
%   escaped. Life is too short to worry about escaping characters, so we just use a file.
tmp_name = tempname();
tmp_fid = fopen(tmp_name, 'w');
tmp_cleanup = onCleanup(@()close_and_delete(tmp_fid));
fprintf(tmp_fid, raw_json);
[status, pretty_json] = system(sprintf('jq . --indent 4 --sort-keys %s', tmp_name));
assert(status == 0, 'External call to jq failed with output: %s', pretty_json);

% append to log file
log_fid = fopen(log_file, 'a+');
log_cleanup = onCleanup(@() fclose(log_fid));
fprintf(log_fid, pretty_json);

end


function close_and_delete(handle) 
    % close the (open) file handle and delete the corresponding file, used for cleanup
    % %
    name = fopen(handle);  % oddly, this is how MATLAB inspects open files
    fclose(handle);
    delete(name);
end
