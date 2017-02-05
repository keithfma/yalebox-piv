function hash = util_md5_hash(file)
% function hash = util_md5_hash(file)
%
% Compute MD5 hash of the input file, otherwise fail with an error. Simple
% wrapper around the FileExchange submission DataHash by Jan Simon.
%
% Arguments:
%
%   file = String, input file to hash
%   hash = String, MD5 hash of input file
% %

opt = struct('Method', 'MD5', 'Format', 'hex', 'Input', 'file');
hash = DataHash(file, opt);
assert(length(hash) == 32, 'Failed to get MD5 hash');