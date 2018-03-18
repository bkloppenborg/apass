#!/bin/bash

# Filename: python-pipeline.sh
# Purpose: Tests the Python pipeline against a subset of the SRO and APASS data
# Usage: Run the script (no arguments required). It'll stop if there is an error


####
# End-to-end test for SRO data
####
export DATA_DIR=../data/other-test-data/
export SAVE_DIR=../data/test-output/
export CODE_DIR=../python/
export OPTS="-j 8"

# clear out the save directory
rm ${SAVE_DIR}/*

# instruct bash to stop on the first non-true exit code
set -e -x
START=$(date +%s)

# run all steps in the SRO pipeline
python ${CODE_DIR}/make_zones.py ${SAVE_DIR}

# test import of malformed FRED file
pdb ${CODE_DIR}/fred_to_zone.py ${OPTS} --debug ${DATA_DIR}/malformed.fred ${SAVE_DIR}

#python ${CODE_DIR}/zone_to_rects.py ${OPTS} ${SAVE_DIR}/*.fredbin
#python ${CODE_DIR}/fix_zone_overlaps.py ${OPTS} ${SAVE_DIR}
#python ${CODE_DIR}/sro_rect_to_dat.py ${OPTS} ${SAVE_DIR}/*-container.fredbin
#python ${CODE_DIR}/sro_merge_dat.py ${SAVE_DIR}/z*.dat

# print out timing statistics
END=$(date +%s)
DIFF=$(( $END - $START ))
echo "SRO pipeline took $DIFF seconds"
