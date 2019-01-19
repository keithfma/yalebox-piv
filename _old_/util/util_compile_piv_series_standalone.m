% Compile piv_series_standalone() wrapper function as standalone executable,
% suitable for use with the SCC batch system.

% generate version-specific output directory
hash = util_git_hash();
output_dir = fullfile('.', 'standalone', hash);

% create / clobber output directory
mkdir(output_dir);

% compile to standalone executable in specified output directory
start_dir = pwd;
cd(output_dir);
mcc('-m', '-R', '-nodisplay', 'piv_series_standalone.m');
cd(start_dir);
