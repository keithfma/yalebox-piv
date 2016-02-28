function hash = util_git_hash()
% function hash = util_git_hash()
%
% Fetch the revision number of the Git repository this file belongs to,
% otherwise fail with error.
% %

git_dir = fileparts(mfilename('fullpath'));
git_cmd = sprintf('git --git-dir %s/.git rev-parse HEAD', git_dir);

[stat, out] = system(git_cmd);
assert(stat == 0, 'Failed to find git revision number');

out_parts = strsplit(out, '\n');
hash = out_parts{1};
assert(length(hash) == 40, 'Failed to extract hash from stdout');