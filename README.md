# Yalebox-PIV

Tools for sandbox image series analysis developed for the Yale University
Department of Geology &amp; Geophysics "Yalebox". Includes:

- `prep`: image pre-processing, e.g.,rectify, crop, mask, etc.
- `piv`: particle image velocimetry
- `post`: displacement vector post-processing, e.g., velocity and strain
- `movie`: image, velocity, strain, etc

Processing steps are dependent, a typical workflow is:

 `prep` &#9658; `piv` &#9658; `post`

_TODO: replace the above with a simple DAG including movies_

For all the details, see *Ma et al 2021 (in prep)*.

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


## Feature flags

During development, it is useful to toggle between different code pathways without either modifiying the
code directly or passing flag parameters all the way down from the top-level interface. We use
["feature flags"](https://www.martinfowler.com/articles/feature-toggles.html) to allow developers to 
experiment with different methods before making a final choice. 

Feature flags are environment variables the user can set to selectively enable/disable features at runtime. 

For this project *all feature flags should be considered temporary*, we will ultimately either choose a single
code path, or "promote" the flag to a normal function argument. If you are adding a feature flag, please document
both here and in the function where it is used, and use the helper function in `subroutines/util/feature_flag.m`
for access and default values.

Current feature flags are:

| Flag name | Description |
| --------- | ----------- |
| YALEBOX\_PIV\_REGISTER\_METHOD | Select which image registration method to use in `piv_register.m` |


## Tests

This package includes unit tests for some core features, but test coverage is by no means complete.

To run the test suite, you can just execute the `yalebox_run_tests.m` script in the root directory.

When writing new tests, please use the [class-based unit tests](https://www.mathworks.com/help/matlab/class-based-unit-tests.html?)
and make sure your new tests are in the set of directories that `yalebox_run_tests.m` searches to 
find tests to run.

