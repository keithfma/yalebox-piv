#!/bin/bash
#
# Utility script for extracting a select number of time steps from a netCDF
# file and saving the results as a new netCDF. This script uses the NCO tool
# NCKS, which does all of the heavy lifting. All this script does is to
# automate the syntax.
# #

# display help message
usage () {
echo "Usage:"
echo "$(basename $0) [input file] [num steps] [variables] [output file]"
echo ""
echo "	input file = path to input netCDF"
echo "	num steps = number of timesteps to extract, these will be equally"
echo "		spaced and span the input dataset as nearly as possible"
echo "	variables = names of variables to extract, as a comma-separated list"
echo "	output file = path to output netCDF"
echo ""
}

# check for sane inputs
if [ $# -ne 4 ]; then
	echo "Incorrect number of arguments"
	usage
fi
input_file=$1
variable_list=$3
output_file=$4


# subset using NCO tool NCKS
ncks -v $variable_list $input_file $output_file 
