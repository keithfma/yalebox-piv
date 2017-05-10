% Compile piv_series_standalone() wrapper function as standalone executable,
% suitable for use with the SCC batch system.

% generate version-specific output directory
os = computer();
release = version('-release');
hash = util_git_hash();
output_dir = fullfile('.', 'standalone', [hash(1:10), '_', os, '_', release]);

% create / clobber output directory
if exist(output_dir, 'dir') == 7
    rmdir(output_dir, 's');
end
mkdir(output_dir);

% compile to standalone executable in specified output directory
start_dir = pwd;
cd(output_dir);
mcc('-m', '-R', '-nodisplay', 'piv_series_standalone.m');
cd(start_dir);
