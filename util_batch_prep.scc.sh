#!/bin/bash -l
#
#  Example submission script for image series pre-processing
#
#$ -N yalebox_prep
#$ -P glaciermod
#$ -j y
#$ -l h_rt=24:00:00

# define parameters
yalebox_path=/projectnb/glaciermod/yalebox-piv
param_file=

# setup runtime environment
module purge
module load matlab/2015b
export MATLABPATH=$yalebox_path:$MATLABPATH

# run
matlab -singleCompThread -nodisplay << EOF
load $param_file
prep_series(output_file, image_path, image_names, ctrl_xw, ctrl_yw, ctrl_xp, ctrl_yp, crop_xw, crop_yw, entropy_len, num_cluster, cluster_center, eql_len, xw, yw, mask_manual);
exit
EOF
