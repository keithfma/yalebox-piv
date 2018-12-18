function [] = test_run(param_file, num_images)
% Run prep for several images in the series and preview results from the
% output netCDF file
%
% Arguments:
%   param_file: string, path to parameter definition JSON file
%   num_images: number of images to include in the test run
% %

param = loadjson(param_file, 'SimplifyCell', 1);

% edit file list

% save new param file

% call prep series with new param file [ REQUIRES REFACTOR FOR PREP_SERIES ]

% call view function to preview results [ REQUIRES NEW VIEW FUNCTION ]