# Yalebox-PIV

Tools for sandbox image series analysis developed for the Yale University
Department of Geology &amp; Geophysics "Yalebox". Includes:

- `prep`: image pre-processing, e.g.,rectify, crop, mask, etc.
- `piv`: particle image velocimetry
- `post`: displacement vector post-processing, e.g., velocity and strain
- `movie`: image, velocity, strain, etc

Processing steps are dependent, a typical workflow is:

 `prep` &#9658; `piv` &#9658; `post`

_TODO: replace the above with a simple DAG_

For all the details, see *Ma et al 2019 (in prep)*.

_TODO: add proper reference when ready_

This package is available under the MIT license. We ask only that you cite the
reference above if you use this in your research. 

## Usage: prep

+ _in_: Experiment images, world-coordinate image, and a bunch of parameters as
  a single JSON file.
+ _out_: Single MAT file with data arrays ready for PIV analysis, self-documenting.

1. Run `prep_template` to generate an initial parameter file. The file is JSON
   and includes help text defining each of the variables.
1. Open the `prep_param` script and run it cell-by-cell. Each cell defines
   and/or test parameters for one pre-processing step. Edit the parameter file
   as needed. We know this is wierd, but it works well. 
1. Run `prep` to execute the pre-processing run you have prepared. This will
   take a while.

## Usage: piv 

+ _in_: output file from prep step, and a bunch of parameters as a single JSON file.
+ _out_: Single MAT file with PIV results, self-documenting.

1. Run `piv_template` to generate an initial parameter file. The file is JSON
   and includes help text defining each of the variables.
1. Open the `piv_param` script and run it cell-by-cell. Each cell defines
   and/or test parameters for one PIV step. Edit the parameter file as
    needed.We know this is wierd, but it works well. 
1. Run `piv` to execute the PIV run you have prepared. This will take a while.

## Usage: post

+ _in_: output file from piv step, and a bunch of parameters as a single JSON file.
+ _out_: Single MAT file with post-processing results, self-documenting.

1. Run `post_template` to generate an initial parameter file. The file is JSON
   and includes help text defining each of the variables.
1. Open the `post_param` script and run it cell-by-cell. Each cell defines
   and/or test parameters for one post-processing step. We know this is wierd,
   but it works well. 
1. Run `post` to execute the post-processing run you have prepared. This will
   take a while.

## Usage: movie

_TODO: movie code needs to be updated_
