Here we test the ability of our PIV routine to resolve offsets in using images
padded with Mark's boundary pattern padding technique. Mark generated the input
image file, which is slightly out of sync with the current file format. The
tasks to complete are:

- [X] Write a script for converting Mark's prepped image file to the current spec
- [X] Run PIV on the converted file
- [ ] Generate figures exploring the results
- [ ] (Optionally) Generate a video from the figures 

Log:

### 2020.10.30

+ wrote `run_piv.m` script to run PIV cases
+ some bugfixes to get the PIV to run, namely (1) new logic for centering the sample grid (the old logic did 
  not make sense to me anymore), and (2) fix the units on the patch image file coordinate vectors.
+ ran a few cases for the PIV, the first (`initial-test`) with just some quick parameters to see if it would run at all,
  and the second (`better-params`) with more reasonable parameters that I would have expected to work. 
+ results are bad outside the wedge we drop basically everything either because we fail to find a peak (common) or
  the validation routine throws the pixel out (more common). Since we drop everything outside the wedge, what is left
  is extrapolated garbage. Saved a figure or two to share with Mark.

### 2020.10.29

+ copied reference image file from `Dropbox/yalebox-piv/experiment/fill-initial-guess/PREP-default-IMAGE.mat` to `images_ref.mat`
+ wrote `patch_image_file.m` to convert source image file to the latest format, and used it to create the `images.mat` file


### 2020.10.28

+ Copied original image file from `Dropbox/Sandbox\ Modeling/RUNS_Phase04\ Faultless\ Wedges/FAULT_SS_01/fault_ss_01_sideN_image.mat` to 
  `images_src.mat`

