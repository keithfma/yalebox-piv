function param_to_json(filename, parameters, clobber)
% function write_params_to_log(log_file, case_name, parameters)
%
% Record case name and input parameters to log file, formatted as concatenated JSON objects.
%
% Makes external system call to `jq` CLI program, which must be installed.
% 
% Arguments:
%   filename: .json file to write serialized parameters to
%   parameters: a yalebox-piv parameters struct (prep, piv, etc all should work)
%   clobber: set True to allow overwriting contents of existing file at 'filename', default
%       is False.
% % 

if nargin < 3; clobber = False; end

% safety checks
assert(endsWith(filename, '.json'), 'Expect filename to have extension .json');
if ~clobber
    assert(isfile(filename) == 0, 'File exists and overwriting is disabled, abort!');
end

% use built-in tools to serialize as JSON
raw_json = jsonencode(parameters, 'ConvertInfAndNaN', true);

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

% write to output file
out_fid = fopen(filename, 'w');
out_cleanup = onCleanup(@() fclose(out_fid));

fprintf(out_fid, pretty_json);

end


function close_and_delete(handle) 
    % close the (open) file handle and delete the corresponding file, used for cleanup
    % %
    name = fopen(handle);  % oddly, this is how MATLAB inspects open files
    fclose(handle);
    delete(name);
end
