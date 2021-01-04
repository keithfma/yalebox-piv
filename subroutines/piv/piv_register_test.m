% unit tests for piv_register()
% %

% TODO: parameterize over valid methods (temporary)

% note: input data used by tests is stored in the file test_data.mat, which contains:
%   + rgb_image: RGB image of sand cropped from a yalebox experiment (fault_ss_01_siden_150.jpg)


classdef piv_register_test < matlab.unittest.TestCase
    
    properties (Constant)
        % TODO: document what this section is for
        % %
        
        % read-only matfile object containing test data
        mat = matfile(...
            fullfile(fileparts(mfilename('fullpath')), 'test_data.mat'), ...
            'Writable', false);
        % normalized source image used to create reference and sample images
        source_image = piv_register_test.load_source_image();
        % row and column coordinate vectors for source image
        source_rows = 1:size(piv_register_test.source_image, 1);
        source_cols = 1:size(piv_register_test.source_image, 2);
    end
    
    properties (TestParameter)
        % TODO: document what this section is for
        % %        
        
        % TODO: describe
        reference_trim = struct(...
            'square', [0, 0], ...
            'tall_odd', [0, 20], ...
            'wide_even', [20, 0]);
        % location of the sample window in the source image, as a [row, column] offset relative to
        %   the center of the source image
        sample_location = struct(...
            'center', [0, 0], ...
            'not_center', [-31.2, 33]);
        % dimensions of sample window
        sample_dims = struct(...
            'square_even', [30, 30], ...
            'square_odd', [31, 31], ...
            'rect_mixed', [30, 41]);
        % TODO: describe
        sample_translation = struct(...
            'none', [0, 0], ...
            'small_int_pos_pos', [10, 12], ... 
            'small_int_pos_neg', [10, -12], ...
            'small_float_pos_pos', [10.1, 11.9], ...
            'small_float_pos_neg', [10.1, -11.9], ...
            'large_float_pos_neg', [70.3, -79.8]);
    end
    
    methods (Static)
        % TODO: document what this section is for
        % %
            
        function img = load_source_image()
            % return normalized image for registration tests
            mat = piv_register_test.mat;
            gray = prep_grayscale(mat.rgb_image);
            img = prep_equalize(gray, true(size(gray)), 31);
        end
        
    end
    
    methods
        % TODO: document what this section is for
        % %        
        
        function [image, rows, cols] = create_reference_image(self, trim)
            % return reference image and its coordinate vectors for specified test case
            % %
            
            validateattributes(trim, {'numeric'}, {'vector', 'integer', 'numel', 2});
            
            self.assertTrue(size(self.source_image, 1) == size(self.source_image, 2), ...
                'expect source image to be square');
            
            image = self.source_image(1:end-trim(1), 1:end-trim(2));
            rows = self.source_rows(1:end-trim(1));
            cols = self.source_cols(1:end-trim(2));
            
            self.assertTrue(all(size(image) == (size(self.source_image) - trim)), 'wrong image size')
            self.assertTrue(length(rows) == size(self.source_image, 1) - trim(1), 'wrong row coord size');
            self.assertTrue(length(cols) == size(self.source_image, 2) - trim(2), 'wrong col coord size');
        end
        
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
            center_idx = (1 + size(self.source_image))/2 + location;  
            row_idx = row_idx + center_idx(1);
            col_idx = col_idx + center_idx(2);
            % ...check results satisfy assumptions
            self.assertEqual(mean(row_idx), center_idx(1), 'AbsTol', 0.51, 'rows not centered');
            self.assertEqual(mean(col_idx), center_idx(2), 'AbsTol', 0.51, 'columns not centered');
            self.assertLength(row_idx, dimensions(1), 'wrong row index length');
            self.assertLength(col_idx, dimensions(2), 'wrong column_index length');
            % note: we check that all indices are within the reference image after interpolation
            
            % resample image and coordinates from reference           
            rows = interp1(1:size(self.source_image, 1), self.source_rows, row_idx, 'linear');
            cols = interp1(1:size(self.source_image, 2), self.source_cols, col_idx, 'linear');
            image = interp2(...
                self.source_cols, ...
                self.source_rows, ...
                self.source_image, ...
                cols, ...
                rows', ...
                'cubic');
            % ...check results satisfy assumptions
            self.assertTrue(all(size(image) == dimensions), 'incorrect dimensions');
            
            % % DEBUG: plot reference image and resampled image on the same axes with markers for 
            % %   the corners of the resampled image
            % figure
            %imagesc(self.source_cols([1, end]), self.source_rows([1, end]), self.source_image);
            % hold on
            % imagesc(cols([1, end]), rows([1, end]), image);
            % colormap('gray');
            % plot(...
            %    [cols(1), cols(1), cols(end), cols(end)], ...
            %     [rows(1), rows(end), rows(end), rows(1)], ...
            %    'Color', 'r', 'Marker', '.', 'MarkerSize', 24, 'LineStyle', 'None');
            % keyboard
            % /DEBUG
            
            % translate image coordinates such that it would have to be translated by 'translation'
            %   to correctly register against the reference image
            rows = rows - translation(1);
            cols = cols - translation(2);
            
        end
        
    end
       
    
    methods (Test)
        % TODO: document what this section is for
        % %        
        
        function test_correct_displacement(self, reference_trim, sample_location, sample_dims, sample_translation)
            % Check if we can recover known displacement with piv_register
            % TODO: define arguments when the dust has settled
            % %
            
            % constants
            ABSOLUTE_TOLERANCE = 0.1;  % pixels
            MIN_OVERLAP = 0.8*prod(sample_dims);
            % /DEBUG
            
            [reference_image, reference_rows, reference_cols] = self.create_reference_image(...
                reference_trim);
            
            [sample_image, sample_rows, sample_cols] = self.create_sample_image(...
                sample_location, sample_dims, sample_translation);
            
            sample_origin = [sample_rows(1), sample_cols(1)];
            sample_mask = true(size(sample_image));
            
            reference_origin = [reference_rows(1), reference_cols(1)];
            reference_mask = true(size(reference_image));
                
            [delta_row, delta_col, quality] = piv_register(...
                sample_image, sample_mask, sample_origin, ...
                reference_image, reference_mask, reference_origin, ...
                MIN_OVERLAP); 
            
            self.verifyEqual(quality, Quality.Valid);
            self.verifyEqual(delta_row, sample_translation(1), 'AbsTol', ABSOLUTE_TOLERANCE);
            self.verifyEqual(delta_col, sample_translation(2), 'AbsTol', ABSOLUTE_TOLERANCE);
            
        end
        
    end
    
end
    