#!/bin/bash
#
#  Example submission script for the BU SCC cluster scheduler
#
#$ -N fault_ss_01_sidef.piv
#$ -P glaciermod
#$ -j y
#$ -pe omp 12

# define parameters
yalebox_path=/projectnb/glaciermod/yalebox-piv
data_top_dir=/projectnb/glaciermod/yalebox-exp-fault/data/fault_01_ss
data_run_name=fault_ss_01_sidef
param_file=$data_top_dir/piv/$data_run_name.piv_param.mat
input_file=$data_top_dir/piv/$data_run_name.input.nc
output_file=$data_top_dir/piv/$data_run_name.displ.nc

# setup runtime environment
module load matlab/2015b
export MATLABPATH=$yalebox_path:$MATLABPATH

# run
matlab -nodisplay << EOF
load $param_file
piv_series('$output_file', '$input_file', samplen, sampspc, intrlen, npass, valid_max, valid_eps, lowess_span_pts, spline_tension, min_frac_data, min_frac_overlap, low_res_spc, verbose)
exit
EOF
