function version_string = get_version()
% ver_str = function get_version()
%
% Return Yalebox version string, read from VERSION.txt
% %

% get version file path, taking care in case called from a weird place
[base_path, ~, ~] = fileparts(mfilename('fullpath'));
version_file = fullfile(base_path, '..', '..',  'VERSION');

% read version string
fid = fopen(version_file);
version_string = fgetl(fid);
fclose(fid);