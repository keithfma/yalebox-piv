function pp = piv_default_param(filename)
% function pp = piv_default_param(filename)
%
% Save PIV parameter file to specified path
% %

update_path('util');

[~, ~, ext] = fileparts(filename);
assert(strcmp(ext, '.json'), 'Filename must have extension .json');

pp.step_range.help = "2-element vector, [initial, final] steps of input images series to include in PIV analysis. Either element can be set to NaN to use the full range.";
pp.step_range.value = [NaN, NaN];

pp.gap.help = "Scalar, integer, gap between analyzed images in steps, for example, for an initial image at step 3, setting gap -> 1 would use a final image at step 4 and yield PIV results at step 3.5, or alternatively, setting gap -> 2 would use a final image at step 5 and yield PIV results at step 4.";
pp.gap.value = 1;

pp.samp_len.help = "Vector, length == number of grid resolutions, integer, side length of the square sample window, [pixels]";
pp.samp_len.value = [30, 30];

pp.samp_spc.help = "Scalar, integer, spacing between adjacent sample points in the (square) sample grid, [pixels].";
pp.samp_spc.value = 30;

pp.intr_len.help = "Vector, length == number of grid resolutions, integer, side length of the square interrogation window";            
pp.intr_len.value = [120, 40];
        
pp.num_pass.help = "Vector, length == number of grid resolutions, integer, number of image deformation passes";
pp.num_pass.value = [1, 2];
       
pp.valid_radius.help = "Scalar, radius around each sample point to include in vector validation, recommended to use ~4*samp_spc, [pixel] units"; 
pp.valid_radius.value = 120;
            
pp.valid_max.help = "Scalar, double, maximum value for the normalized residual in the vector validation function, above which a vector is flagged as invalid. Ref [3] reccomends a value of 2.";
pp.valid_max.value = 2;
            
pp.valid_eps.help = "Scalar, double, minimum value of the normalization factor in the vector validation function. Ref [3] reccomends a value of 0.1.";
pp.valid_eps.value = 0.1;
            
pp.min_frac_data.help = "Scalar, minimum fraction of the sample window that must contain data (e.g. sand) for the point to be included in the ROI for PIV analysis";
pp.min_frac_data.value = 0.8;
            
pp.min_frac_overlap.help = "Scalar, minimum fraction of the sample window data that must overlap the interrogation window data for a point in the cross-correlation to be valid"; 
pp.min_frac_overlap.value = 0.5;

save_param(pp, filename);
