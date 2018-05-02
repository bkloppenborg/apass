#!/bin/bash

# Filename: run-pipeline.sh
# Purpose: Runs the correlation and averaging stage of the APASS pipeline.
# Usage: Run the script (no arguments required). It'll stop if there is an error

####
# Directory Configuration
####
#export DATA_DIR=/mnt/astro/astro/home/welch/surveys/APASS/hphot_and_hepp/
export DATA_DIR=/mnt/astro/astro/home/welch/surveys/APASS/hphot_and_hepp/fixedk2_dtfred/
export CACHE_DIR=/scratch/apass-dr10/
export CODE_DIR=./python/
export SAVE_DIR=/home/bkloppenborg-gmail_com/apass-dr10/
export OPTS="-j 45" # Astro has 24 cores, 48 threads. Leave one core free.

# clear out the cache directory
mkdir -p ${CACHE_DIR}
rm ${CACHE_DIR}/*

# clear out the save directory
mkdir -p ${SAVE_DIR}
rm ${SAVE_DIR}/*

# instruct bash to stop on the first non-true exit code
set -e -x
START=$(date +%s)

# run all steps in the pipeline.
# sync files to the save directory
python ${CODE_DIR}/make_zones.py ${CACHE_DIR}
python ${CODE_DIR}/fred_to_zone.py ${OPTS} ${DATA_DIR}/*.fred ${CACHE_DIR}
python ${CODE_DIR}/zone_to_rects.py ${OPTS} ${CACHE_DIR}/*.fredbin
python ${CODE_DIR}/fix_zone_overlaps.py ${OPTS} ${CACHE_DIR}
python ${CODE_DIR}/rect_to_dat.py ${OPTS} apass ${CACHE_DIR}/*-container.fredbin

# print out timing statistics
END=$(date +%s)
DIFF=$(( $END - $START ))
echo "APASS pipeline took $DIFF seconds"

# Move the data to the save directory, clear out the cache directory
echo "Moving data to save directory."
START=$(date +%s)
rsync -avh ${CODE_DIR} ${SAVE_DIR}
DIFF=$(( $END - $START ))
echo "Synchronization took ${DIFF} seconds"
rm -rf ${CACHE_DIR}
