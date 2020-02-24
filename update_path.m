function [] = update_path(varargin)
% [] = function update_path(varargin)
% 
% Load specified dependencies to the matlab path
%
% Arguments: string(s), specific named packages to add to path, invalid
%   names raise an error
%
% Returns: Nothing
% % 

[base_path, ~, ~] = fileparts(mfilename('fullpath'));
dependencies = fullfile(base_path, 'dependencies');
subroutines = fullfile(base_path, 'subroutines');

for ii = 1:length(varargin)
    
    switch varargin{ii}
        
        % core
        case 'main'
            addpath(base_path);
        
        % third-party libraries
        case 'datahash'
            addpath(fullfile(dependencies, 'DataHash_20181113'));
        case 'normxcorr2_masked'
            addpath(fullfile(dependencies, 'MaskedFFTRegistrationCode'));
        case 'spline'
            addpath(fullfile(dependencies, 'spline2d'));
        case 'deriv'
            addpath(fullfile(dependencies, 'derivatives'));
        case 'akde'
            addpath(fullfile(dependencies, 'akde'));
        case 'inpaint_nans'
            addpath(fullfile(dependencies, 'inpaint_nans'));
        case 'export_fig'
            addpath(fullfile(dependencies, 'export_fig'));
        case 'allcomb'
            addpath(fullfile(dependencies, 'allcomb'));
        case 'newmatic'
            addpath(fullfile(dependencies, 'newmatic'));
        
        % subroutines
        case 'prep'
            addpath(fullfile(subroutines, 'prep'));
        case 'piv'
            addpath(fullfile(subroutines, 'piv'));
        case 'post'
            addpath(fullfile(subroutines, 'post'));
        case 'util'
            addpath(fullfile(subroutines, 'util'));
        case 'plot'
            addpath(fullfile(subroutines, 'plot'));
        case 'movie'
            addpath(fullfile(subroutines, 'movie'));
        
        % bad name, fail
        otherwise
            error('Invalid dependency name: "%s"', varargin{ii});
        
    end
end

end
