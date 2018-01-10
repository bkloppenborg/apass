#!/bin/bash

# Tests the pipeline against a subset of SRO and APASS data


# End-to-end test for SRO data
export DATA_DIR=/home/data/sro-test-data/
export SAVE_DIR=/home/data/sro-test/

python make_zones.py ${SAVE_DIR}
python fred_to_zone.py -j 8 ${DATA_DIR}/*.fred ${SAVE_DIR}
python zone_to_rects.py -j 8 ${SAVE_DIR}/*.fredbin
python fix_zone_overlaps.py -j 8 ${SAVE_DIR}
python sro_rect_to_dat.py -j 8 ${SAVE_DIR}/*-container.fredbin
python sro_merge_dat.py ${SAVE_DIR}/z*.dat
