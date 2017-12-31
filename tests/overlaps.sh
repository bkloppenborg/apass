#!/bin/bash

# Filename: overlaps.sh
# Purpose: Tests the Python pipeline against a fake data that is crafted to
#          overlap within zones and between zones.
# Usage: Run the script (no arguments required). It'll stop if there is an error


####
# End-to-end test to check for overlapping data points.
####
export DATA_DIR=../data/overlap-test-data/
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
python ${CODE_DIR}/fred_to_zone.py ${OPTS} ${DATA_DIR}/*.fred ${SAVE_DIR}
python ${CODE_DIR}/zone_to_rects.py ${OPTS} ${SAVE_DIR}/*.fredbin
python ${CODE_DIR}/fix_zone_overlaps.py ${OPTS} ${SAVE_DIR}
python ${CODE_DIR}/sro_rect_to_dat.py ${OPTS} ${SAVE_DIR}/*-container.fredbin

# print out timing statistics
END=$(date +%s)
DIFF=$(( $END - $START ))
echo "Pipeline took $DIFF seconds"

# plot the four zones
python ${CODE_DIR}/plot_zone.py ${SAVE_DIR}/z03351.fredbin --show-containers &
python ${CODE_DIR}/plot_zone.py ${SAVE_DIR}/z03394.fredbin --show-containers &
python ${CODE_DIR}/plot_zone.py ${SAVE_DIR}/z01985.fredbin --show-containers &
python ${CODE_DIR}/plot_zone.py ${SAVE_DIR}/z02028.fredbin --show-containers &
