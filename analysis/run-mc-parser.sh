#!/bin/sh

# shell script to run mc_parse.C for use with submit-mc-parser-jobs.sh

## Usage
#./run-mc-parser.sh <config-file-name> 

file=$1

myfile='"'$file'"'

cd /w/halla-scshelf2102/sbs/ktevans/KateJackSBSAnalysis/analysis

root -l -b -q 'mc_parse.C('$myfile')'
