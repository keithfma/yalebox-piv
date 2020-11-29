Here we test the ability of our PIV routine to resolve offsets in using images
padded with Mark's boundary pattern padding technique. Mark generated the input
image file, which is slightly out of sync with the current file format. The
tasks to complete are:

- [X] Write a script for converting Mark's prepped image file to the current spec
- [ ] Run PIV on the converted file
- [ ] Generate figures exploring the results
- [ ] (Optionally) Generate a video from the figures 

Log:

### 2020.10.29

+ copied reference image file from `Dropbox/yalebox-piv/experiment/fill-initial-guess/PREP-default-IMAGE.mat` to `images_ref.mat`
+ wrote `patch_image_file.m` to convert source image file to the latest format, and used it to create the `images.mat` file


### 2020.10.28

+ Copied original image file from `Dropbox/Sandbox\ Modeling/RUNS_Phase04\ Faultless\ Wedges/FAULT_SS_01/fault_ss_01_sideN_image.mat` to 
  `images_src.mat`

