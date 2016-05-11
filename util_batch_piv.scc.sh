#!/bin/bash
#
#  Example submission script for the BU SCC cluster scheduler
#
#$ -N K23_side_piv
#$ -P glaciermod
#$ -j y
#$ -pe omp 12
#$ -l h_rt=24:00:00

# define parameters
yalebox_path=/projectnb/glaciermod/yalebox-piv
param_file=/projectnb/glaciermod/yalebox-exp-erosion/data/K23_side_piv_param.mat
input_file=/projectnb/glaciermod/yalebox-exp-erosion/data/K23_side.image.nc
output_file=/projectnb/glaciermod/yalebox-exp-erosion/data/K23_side.piv.nc

# setup runtime environment
module load matlab/2015b
export MATLABPATH=$yalebox_path:$MATLABPATH

# run
matlab -nodisplay << EOF
load $param_file
util_setup_env();
piv_series('$output_file', '$input_file', samplen, sampspc, intrlen, npass, valid_max, valid_eps, lowess_span_pts, spline_tension, min_frac_data, min_frac_overlap, true);
exit
EOF
