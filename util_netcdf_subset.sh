#!/bin/bash
#
# Utility script for extracting a select number of time steps from a netCDF
# file and saving the results as a new netCDF. This script uses the NCO tool
# NCKS, which does all of the heavy lifting. All this script does is to
# automate the syntax.
# #

usage () {
echo "Usage:"
echo "$(basename $0) [input file] [num steps] [variables]"
echo ""
echo "	input file = path to input netcdf"
echo "	num steps = number of timesteps to extract, these will be equally"
echo "		spaced and span the input dataset as nearly as possible"
echo "	variables = names of variables to extract, as a comma-separated list"
echo ""
}

if [ $# -ne 3 ]; then
	echo "Incorrect number of arguments"
	usage
fi
