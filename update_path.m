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


function [] = add_dep(dep_path)
    addpath(dep_path);
    fprintf('Added %s to MATLAB path\n', dep_path);
end

for ii = 1:length(varargin)
    
    switch varargin{ii}
        
        % third-party libraries
        case 'jsonlab'
            add_dep(fullfile(dependencies, 'jsonlab-1.8'));
        case 'xcorr'
            add_dep(fullfile(dependencies, 'MaskedFFTRegistrationCode'));
    
        % subroutines
        case 'prep'
            add_dep(fullfile(subroutines, 'prep'));
        case 'piv'
            add_dep(fullfile(subroutines, 'piv'));
        case 'util'
            add_dep(fullfile(subroutines, 'util'));
        
        % bad name, fail
        otherwise
            error('Invalid dependency name: "%s"', dep_name);
        
    end
end

end