#!/bin/bash
#
#  Submission script for experiment fault_ss_02_siden PIV
#
#$ -N fault_ss_02_siden_piv
#$ -P glaciermod
#$ -j y
#$ -pe omp 8 
#$ -l h_rt=72:00:00

# define parameters
yalebox_path=/projectnb/glaciermod/yalebox-piv
param_file=/projectnb/glaciermod/yalebox-exp-fault/data/fault_ss_02_siden_piv_param.mat
input_file=/projectnb/glaciermod/yalebox-exp-fault/data/fault_ss_02.siden_image.nc
output_file=/projectnb/glaciermod/yalebox-exp-fault/data/fault_ss_02_siden_piv.nc

# setup runtime environment
module load matlab/2015b
export MATLABPATH=$yalebox_path:$MATLABPATH

# run
matlab -nodisplay << EOF
load $param_file
util_setup_env();
maxNumCompThreads(getenv('NSLOTS'));
fprintf('Maximum number of threads is set to %d\n', maxNumCompThreads);
piv_series('$output_file', '$input_file', samplen, sampspc, intrlen, npass, valid_max, valid_eps, lowess_span_pts, spline_tension, min_frac_data, min_frac_overlap, true);
exit
EOF
