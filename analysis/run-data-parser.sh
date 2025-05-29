#!/bin/sh

# shell script to run data_elastic_parse.C for use with submit-data-parser-jobs.sh

## Usage
#./run-data-parser.sh <config-file-name> 

file=$1

myfile='"'$file'"'

cd /w/halla-scshelf2102/sbs/ewertz/nTPE_analysis/analysis/gmn_ntpe

root -l -b -q 'data_elastic_parse.C('$myfile')'
