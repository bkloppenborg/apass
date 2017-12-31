SRO Test Data
=====

This directory contains a few files which serve as test data for the SRO portion
of the APASS pipeline. The data consist of the first 100 lines from each SRO FRED
file as of 2017-12-31. The files were generated as follows:

    cd sro-data
    mkdir sro-test-data
    for f in *.fred; do head -n 100 $f > ./sro-test-data/$f; done

Data generated in this fashion are likely sufficient to test most of the pipeline
**except** the `fix_zone_overlap.py` script as it is unlikely that data from two
stars will be split across two different zones. As such, it may be necessary to
manually create this scenario using fictitious data.
