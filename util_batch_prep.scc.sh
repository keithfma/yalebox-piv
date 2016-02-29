#!/bin/bash
#
#  Example submission script for image series pre-processing
#
#$ -N fault_ss_01_sidef.prep
#$ -P glaciermod
#$ -j y
#$ -l h_rt=24:00:00

# define parameters
yalebox_path=/projectnb/glaciermod/yalebox-piv
data_top_dir=/projectnb/glaciermod/yalebox-exp-fault/data/fault_ss_01
data_run_name=fault_ss_01_sidef
param_file=$data_top_dir/piv/$data_run_name.prep_param.mat
output_file=$data_top_dir/piv/$data_run_name.input.nc
image_path=$data_top_dir/image/crop_sidef
image_wild=${image_path}'/fault_ss_01_sidef_*.png'

# setup runtime environment
module load matlab/2015b
export MATLABPATH=$yalebox_path:$MATLABPATH

# run
matlab -singlecompthread -nodisplay << EOF
load $param_file
d = dir('$image_wild'); image_names = {d.name};
prep('$output_file', '$image_path', image_names, x, y, scale, offset, mask_manual, hue_lim, val_lim, entr_lim, entr_win, morph_open_rad, morph_erode_rad, nwin);
exit
EOF
