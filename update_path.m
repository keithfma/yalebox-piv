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
        
        % third-party libraries
        case 'jsonlab'
            addpath(fullfile(dependencies, 'jsonlab-1.8'));
        case 'datahash'
            addpath(fullfile(dependencies, 'DataHash_20181113'));
        case 'normxcorr2_masked'
            addpath(fullfile(dependencies, 'MaskedFFTRegistrationCode'));
        case 'spline'
            addpath(fullfile(dependencies, 'spline2d'));
        case 'deriv'
            addpath(fullfile(dependencies, 'derivatives'));
    
        % subroutines
        case 'prep'
            addpath(fullfile(subroutines, 'prep'));
        case 'piv'
            addpath(fullfile(subroutines, 'piv'));
        case 'post'
            addpath(fullfile(subroutines, 'post'));
        case 'util'
            addpath(fullfile(subroutines, 'util'));
        
        % bad name, fail
        otherwise
            error('Invalid dependency name: "%s"', dep_name);
        
    end
end

end