function hash = util_md5_hash(file)
% function hash = util_md5_hash(file)
%
% Compute MD5 hash of the input file, otherwise fail with an error. Currently
% uses a system call to an external system tool (md5sum, linux).
%
% Arguments:
%
%   file = String, input file to hash
%   hash = String, MD5 hash of input file
% %

[stat, out] = system(sprintf('md5sum %s', file));
assert(stat == 0, 'Failed to compute MD5 hash');

out_parts = strsplit(out, ' ');
hash = out_parts{1};
assert(length(hash) == 32, 'Failed to extract value of MD5 hash from stdout');