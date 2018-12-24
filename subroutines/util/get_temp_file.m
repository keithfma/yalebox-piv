function filename = get_temp_file(extension)
% function filename = get_tempfile(extension)
%
% Return unique filename in yalebox-piv tmp directory
%
% Arguments:
%   extension: string, filename extension 
% %

[base_path, ~, ~] = fileparts(mfilename('fullpath'));
tmp_dir = fullfile(base_path, '..', '..', 'tmp');
filename = [tempname(tmp_dir), '.', strip(extension, 'left', '.')];
