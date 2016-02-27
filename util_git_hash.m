function hash = get_git_hash()
% Fetch the revision number of the Git repository this file belongs to,
% otherwise fail with error.
% %

git_dir = fileparts(mfilename('fullpath'));
git_cmd = sprintf('git --git-dir %s/.git rev-parse HEAD', git_dir);
[stat, hash] = system(git_cmd);
assert(stat == 0, 'Failed to find git revision number');

end 