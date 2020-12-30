% unit tests for piv_register()
% %

% TODO: parameterize over valid methods (temporary)

% TODO: write a function that uses the sand image to generate sample and interrogation windows
%   with suitable coordinate vectors for each, given, ...
%   + size and position (in normalized coordinates) of the sample window
%   + added translation (applied to coordinates only)

% TODO: test:
%   + integer offsets in x, y and both
%   + non-integer offsets in x, y, and both
%   + sample image located in all six "directions" in the interogation image


% note: input data used by tests is stored in the file test_data.mat, which contains:
%   + rgb_image: RGB image of sand cropped from a yalebox experiment (fault_ss_01_siden_150.jpg)


classdef piv_register_test < matlab.unittest.TestCase
    
    properties (Constant)
        image = normalized_source_image()
    end
    
    methods (Test)
        
        function test_nothing(self)
            % Dummy test to see if we can get tests to run at all
            % TODO: delete me!
            fprintf('Did nothing successfully\n'); 
        end
        
        function test_correct_displacement(self)
            % Check if we can recover known displacement with piv_register
            
            reference_image = self.nimage;
            reference_rows = 1:size(reference_image, 1);
            reference_cols = 1:size(reference_image, 2);
            
            % TODO NEXT: generate sample image and it's coordinates
            keyboard
            
        end
        
    end
    
end


function obj = test_data_mat()
    % return read-only matfile object for the test data MAT-file
    obj = matfile(...
        fullfile(fileparts(mfilename('fullpath')), 'test_data.mat'), ...
        'Writable', false);
end


function norm_image = normalized_source_image()
    % return normalized image for registration tests
    mat = test_data_mat();
    rgb_image = mat.rgb_image;
    gray_image = prep_grayscale(rgb_image);
    norm_image = prep_equalize(gray_image, true(size(gray_image)), 31);
end