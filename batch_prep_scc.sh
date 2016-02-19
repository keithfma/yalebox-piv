#!/bin/bash
#
#  Template submission script for the BU SCC cluster scheduler
#
#$ -N yalebox_prep
#$ -P glaciermod
#$ -j y

# define parameters
yalebox_path=/projectnb/glaciermod/yalebox-piv/
param_file=/projectnb/glaciermod/sandbox-tmp/prep/fault_ss_01_siden_param.mat
input_file_work=$TMPDIR/tmp.nc
input_file_dest=/projectnb/glaciermod/sandbox-tmp/prep/fault_ss_01_siden_in.nc
image_path=/projectnb/glaciermod/sandbox-tmp/fault_ss_01/siden/
image_wild=$image_path'fault_ss_01_siden_*.png'

# setup runtime environment
module load matlab/2014b
export MATLABPATH=$yalebox_path:$MATLABPATH

# run
matlab -singleCompThread -nodisplay << EOF
load $param_file
d = dir('$image_wild'); image_names = {d.name};
prep('$input_file_work', '$image_path', image_names, ...
  x, y, scale, offset, mask_manual, hue_lim, val_lim, entr_lim, entr_win, ...
  morph_rad, num_tiles);
movefile('$input_file_work', '$input_file_dest');
exit
EOF

