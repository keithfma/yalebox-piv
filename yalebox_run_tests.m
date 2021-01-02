% Run all available unit tests
% %

fprintf('Yalebox-PIV unit tests: \n\n');

update_path('prep', 'piv', 'util', 'normxcorr2_masked');

yalebox_test_results = runtests( ...
    'subroutines/piv', ...   % note: could use a cell array if we need to search multiple paths
    'IncludeSubfolders', true, ...
    'OutputDetail', matlab.unittest.Verbosity.Concise ...
);

fprintf('Results summary:\n');

disp(yalebox_test_results);
