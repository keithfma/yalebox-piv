# Yalebox-PIV

Tools for sandbox image series analysis developed for the Yale University
Department of Geology &amp; Geophysics "Yalebox". Includes:

- Preprocessing: distortion correction, masking, adaptive histogram equalization 
- Particle Image velocimetry (PIV): iterative image-deformation method using masked cross-correlation
- Strain analysis
- Movie composition

## Usage

A typical workflow is to make use of the "helper" scripts `prep_get_param.m`
and `piv_get_param.m` to experiment with free parameters for a particular
experiment image series. Once a satisfactory set of parameters is found, the
`prep_series()` and `piv_series()` functions can be used to process the whole
experiment. Results are saved in as internally documented netCDF files. 

## To Do

+ Experiment: use interpolant rather than LOESS smoother to get regular grid
+ Experiment: remove scattered points that lie outside the future footprint
