#!/bin/bash
#
#  Template submission script for the BU SCC cluster scheduler
#
#$ -N yalebox_prep
#$ -P glaciermod
#$ -j y

# define parameters
yalebox_path=/projectnb/glaciermod/yalebox-piv/
param_file=/projectnb/glaciermod/yalebox-exp-fault/data/fault_ss_01/piv/test/prep_fault_ss_01_sidef.mat
output_file_work=/projectnb/glaciermod/yalebox-exp-fault/data/fault_ss_01/piv/test/tmp.nc
output_file_dest=/projectnb/glaciermod/yalebox-exp-fault/data/fault_ss_01/piv/test/out.nc
image_path=/projectnb/glaciermod/yalebox-exp-fault/data/fault_ss_01/piv/test/
image_wild=$image_path'fault_ss_01_sidef_*.png'

# setup runtime environment
module load matlab/2015b
export MATLABPATH=$yalebox_path:$MATLABPATH

# run
matlab -singlecompthread -nodisplay << EOF
load $param_file
d = dir('$image_wild'); image_names = {d.name};
prep('$output_file_work', '$image_path', image_names, x, y, scale, offset, mask_manual, hue_lim, val_lim, entr_lim, entr_win, morph_open_rad, morph_erode_rad, nwin);
movefile('$output_file_work', '$output_file_dest');
exit
EOF

