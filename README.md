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

1. Boundary treatment for derivatives:

Read up on how spline-in-tension interpolation treats boundaries, and on the
details of the derivative filter. It may make sense to do some extrapolation
to prepare for gradient estimation.

1. Explore using vector median (Liu) for rejection filter

1. Explore finite strain analysis

I strongly suspect that treating the displacements as infinitesimal strain is
not valid. More importantly, I suspect that integrating finite strain over a
few time steps will greatly reduce the noise in the resulting figures.
