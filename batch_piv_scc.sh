#!/bin/bash
#
#  Template submission script for the BU SCC cluster scheduler
#
#$ -N yalebox_piv
#$ -P glaciermod
#$ -j y
#$ -pe omp 12

# define parameters
yalebox_path=/projectnb/glaciermod/yalebox-piv/
param_file=/projectnb/glaciermod/yalebox-exp-fault/data/
input_file=/projectnb/glaciermod/yalebox-exp-fault/data/
output_file=/projectnb/glaciermod/yalebox-exp-fault/data/

# setup runtime environment
module load matlab/2015b
export MATLABPATH=$yalebox_path:$MATLABPATH

# run
matlab -nodisplay << EOF
load $param_file
piv_series($output_file, $input_file, samplen, sampspc, intrlen, npass, valid_max, valid_eps, lowess_span_pts, spline_tension, min_frac_data, min_frac_overlap, low_res_spc, verbose)
exit
EOF
