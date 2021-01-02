% unit tests for piv_register()
% %

% TODO: parameterize over valid methods (temporary)

% TODO: test:
%   + integer offsets in x, y and both
%   + non-integer offsets in x, y, and both
%   + sample image located in all six "directions" in the interogation image


% note: input data used by tests is stored in the file test_data.mat, which contains:
%   + rgb_image: RGB image of sand cropped from a yalebox experiment (fault_ss_01_siden_150.jpg)


classdef piv_register_test < matlab.unittest.TestCase
    
    properties (Constant)
        % read-only matfile object containing test data
        mat = matfile(...
            fullfile(fileparts(mfilename('fullpath')), 'test_data.mat'), ...
            'Writable', false);
        % normalized reference image
        reference_image = piv_register_test.create_reference_image();
        
        % row and column coordinate vectors for reference image
        reference_rows = 1:size(piv_register_test.reference_image, 1);
        reference_cols = 1:size(piv_register_test.reference_image, 2);
    end
    
    methods (Static)
            
        function img = create_reference_image()
            % return normalized image for registration tests
            mat = piv_register_test.mat;
            gray = prep_grayscale(mat.rgb_image);
            img = prep_equalize(gray, true(size(gray)), 31);
        end
        
    end
    
    methods
        
        function [image, rows, cols] = create_sample_image(self, location, dimensions, translation)
            % return sample image and its coordinate vectors for specified test case
            % TODO: define arguments when the dust has settled
            % %
            
            validateattributes(dimensions, {'numeric'}, {'vector', 'numel', 2});
            
            % get row and column indices at which we will resample the sample image from the reference
            % note: these "indices" need not be integers, we will resample the image and its coords
            % ...start with index vectors centered about zero
            row_idx = (1:dimensions(1)) - dimensions(1)/2;
            col_idx = (1:dimensions(2)) - dimensions(2)/2;
            % ...adjust vectors so the new center is the reference center + location offset 
            center_idx = (1 + size(self.reference_image))/2 + location;  
            row_idx = row_idx + center_idx(1);
            col_idx = col_idx + center_idx(2);
            % ...check results satisfy assumptions
            self.assertEqual(mean(row_idx), center_idx(1), 'AbsTol', 0.5, 'rows not centered');
            self.assertEqual(mean(col_idx), center_idx(2), 'AbsTol', 0.5, 'columns not centered');
            self.assertLength(row_idx, dimensions(1), 'wrong row index length');
            self.assertLength(col_idx, dimensions(2), 'wrong column_index length');
            % note: we check that all indices are within the reference image after interpolation
            
            % resample image and coordinates from reference
            rows = interp1(1:size(self.reference_image, 1), self.reference_rows, row_idx, 'linear');
            cols = interp1(1:size(self.reference_image, 2), self.reference_rows, row_idx, 'linear');
            image = interp2(...
                1:size(self.reference_image, 1), ...
                1:size(self.reference_image, 2), ...
                self.reference_image, ...
                row_idx', ...
                col_idx, ...
                'cubic');
            
            % DEBUG: plot reference image and resampled image on the same axes with markers for 
            %   the corners of the resampled image
            figure
            imagesc(self.reference_rows([1, end]), self.reference_cols([1, end]), self.reference_image);
            hold on
            imagesc(rows([1, end]), cols([1, end]), image);
            colormap('gray');
            plot(...
                [rows(1), rows(end), rows(end), rows(1)], ...
                [cols(1), cols(1), cols(end), cols(end)], ...
                'Color', 'r', 'Marker', '.', 'MarkerSize', 24, 'LineStyle', 'None');
            % /DEBUG
            
            % translate image coordinates such that it would have to be translated by 'translation'
            %   to correctly register against the reference image
            rows = rows - translation(1);
            cols = cols - translation(2);
            
        end
        
    end
       
    
    methods (Test)
        
        function test_correct_displacement(self)
            % Check if we can recover known displacement with piv_register
            
            % DEBUG: constant versions of parameters
            samp_location = [50, 50];
            samp_dims = [31, 31];
            samp_translation = [0, 0];
            min_overlap = 0.8*prod(samp_dims);
            % /DEBUG
            
            [sample_image, sample_rows, sample_cols] = self.create_sample_image(...
                samp_location, samp_dims, samp_translation);
            
            sample_origin = [sample_rows(1), sample_cols(1)];
            sample_mask = true(size(sample_image));
            
            reference_origin = [self.reference_rows(1), self.reference_cols(1)];
            reference_mask = true(size(self.reference_image));
            
            
            [delta_row, delta_col, quality] = piv_register(...
                sample_image, sample_mask, sample_origin, ...
                self.reference_image, reference_mask, reference_origin, ...
                min_overlap); 
            
            self.assertEqual(quality, Quality.Valid);
            self.assertEqual(delta_row, samp_translation(1), 'AbsTol', 0.01);
            self.assertEqual(delta_col, samp_translation(2), 'AbsTol', 0.01);
            
        end
        
    end
    
end
    