% unit tests for piv_register()
% %

% TODO: parameterize over valid methods (temporary)


% TODO: make a fixture that loads a few samples of sand image and normalizes them, or, loads
%   a pre-normalized sand image. I prefer the former approach because we stay current.

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
        test_data_mat = matfile(...
            fullfile(fileparts(mfilename('fullpath')), 'test_data.mat'), ...
            'Writable', false);
        rgb_image = piv_register_test.test_data_mat.rgb_image;
        gray_image = prep_grayscale(piv_register_test.rgb_image);
        norm_image = prep_equalize(...
            piv_register_test.gray_image, ...
            true(size(piv_register_test.gray_image)), ...
            31);
    end
    
    methods (Test)
        
        function test_nothing(self)
            % Dummy test to see if we can get tests to run at all
            fprintf('Did nothing successfully\n'); 
        end
        
        function test_correct_displacement(self)
            % Check if we can recover known displacement with piv_register
            
            reference_image = self.norm_image;
            reference_rows = 1:size(reference_image, 1);
            reference_cols = 1:size(reference_image, 2);
            
            % TODO NEXT: generate sample image and it's coordinates
            
        end
        
    end
    
end


function img = load_and_normalize_source_image(f)
    img = 5;
end