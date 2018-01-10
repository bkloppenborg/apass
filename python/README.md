APASS Python Pipeline
=====

The contents of this directory represent the final stages of the APASS pipeline.
Here we take the nightly extracted photometry files (.FRED) and reduce them
to a single average measurement for each star in each filter.

## The big picture

The pipeline has the following stages:

* Creation of zone definition files (`make-zones.py`)
* Import star measurement (.FRED) files and split them into zones
  (`fred-to-zone.py`)
* Establish bounding boxes for each star and populate the boxes with
  photometric data (`zone-to-rects.py`)
* Resolve any issues where star measurements overlap between zones
  (`fix-zone-overlaps.py`)
* Perform statistical analysis on the photometric data to generate average
  measurements for each star (`rect-to-data.py`)
  
Most of the Python scripts take command line arguments which can be listed
by specifying the `-h` flag to each program. Additionally, all scripts except
`fix-zone-overlaps.py` can execute in parallel using the `-j` flag.

There are a few configuration variables that are set in the `apass.py` file
which should be changed prior to execution. See below for further information.

## Details


## Prerequisites

* pip
* numpy
* argparse
* networkx 2.0

These packages can be installed on a Ubuntu system with the following commands:

    sudo apt install python-pip python-numpy python-matplotlib
    pip install numpy argparse matplotlib python-tk
    pip install pydotplus
    pip install 'networkx==1.11'


## Sample reduction run

Below is a worked example of running the data reduction pipeline on a subset of
the APASS DR9 FRED files. Sample output from each programs is listed in addition
to contents found in the output directories.

For the purposes of this sample, the following directories are used:

    APASS Code:           ~/workspace/apass/python 
    FRED directory:       /home/data/sro-test-data
    APASS save directory: /home/data/sro-test
    
Create zones using `make_zones.py`

    ~/workspace/apass/python$ python make_zones.py
    Created 4096 zones

This will result in a `global.json` file in the APASS save directory:

    ~/workspace/apass/python$ ls /home/data/apass-test/
    global.json
   
Next we parse each FRED file and split the data into different zones. To run
this process, we use the `fred-to-zone.py` script with `.fred` files as input:

    ~/workspace/apass/python$ python fred_to_zone.py -j 8 /home/data/sro-test-data/*.fred
    Processing FRED file /home/data/sro-test-data/141111.fred
    Processing FRED file /home/data/sro-test-data/141127.fred
    Processing FRED file /home/data/sro-test-data/141024.fred
    Processing FRED file /home/data/sro-test-data/141108.fred
    Processing FRED file /home/data/sro-test-data/141110.fred
    Processing FRED file /home/data/sro-test-data/141126.fred
    ...
    A list of modified files has been written to /home/data/sro-test/fred-to-zone-modified-files.txt
    
This command will generate a series of `.fredbin` files along with a
`*-contrib.txt` file as seen in the output directory below:

    /home/data/sro-test$ ls | head -n 6
    fred-to-zone-modified-files.txt
    global.json
    z00650-contrib.txt
    z00650.fredbin
    z00652-contrib.txt
    z00652.fredbin
  
The `.fredbin` files are a binary representation of the `.fred` files with a few
additonal columns (for bookkeeping zone, node, and container IDs) appended to
the end. If you wish to see the content of the `.fredbin`, there is a `dump-fredbin.py`
script included with the data reduction pipeline. The `.contrib` files describe
which FRED files contributed to data in that zone.

The `fred-to-zone-modified-files.txt` file contains a list of all fredbin files
modified by the execution of `fred-to-zone.py`. You can use bash to expand the
contents of that file to feed in to the next stage of the pipeline as follows:

    python zone-to-rects.py $(< fred-to-zone-modified-files.txt)

Next we parse the zone data and insert the data into a quadtree data structure,
separating the data into discrete bins called "containers". These containers will
(ideally) store measurements for precisely one star, but due to limited resolution
of the telescope and astrometric errors, it is possible that more than one star's
data will end up in a container. These blends will be addressed in a later step.
This stage of the pipeline is completed using the `zone-to-rects.py` script:

    ~/workspace/apass/python$ python zone_to_rects.py -j 8 /home/data/sro-test/*.fredbin
    Processing '/home/data/sro-test/z00705.fredbin' which has 12 data points 
    Processing '/home/data/sro-test/z00652.fredbin' which has 3595 data points 
    Processing '/home/data/sro-test/z00650.fredbin' which has 1566 data points 
    Processing '/home/data/sro-test/z00727.fredbin' which has 370 data points 
    ...
    A list of modified files has been written to /home/data/sro-test/zone-to-rects-modified-files.txt

The zone-to-rect step of the pipeline is the longest operation of all of the
pipeline stages. You can use the `zone-to-rects-modified-files.txt` to supply
files to the `rect-to-dat` stage of the pipeline using the `$(< filename)`
expansion  operator.

After the zone-to-rects script completes, you will find directory contents 
similar to the following:
 
    /home/data/sro-test$ ls | head -n 11
    fred-to-zone-modified-files.txt
    global.json
    z00650-border-rects.json
    z00650-container.fredbin
    z00650-contrib.txt
    z00650.fredbin
    z00650-zone.json
    z00652-border-rects.json
    z00652-container.fredbin
    z00652-contrib.txt
    z00652.fredbin

The `*-container.fredbin` files contain the same data as the previous `.fredbin`
files, but are now sorted/grouped by star and contain non-zero entries for the
node and container IDs. 
The script also generates a series of `*-border-rects.json` files which 
describe any containers whose data might span more than one zone.
    
Next we check for any overlaps between zones using `fix-zone-overlaps.py`. At
present this step does not have a parallel execution option. It automatically
parses data found in the APASS save directory, so simply run the script as
follows:

    ~/workspace/apass/python fix-zone-overlaps.py -j 8 
    Processing polar zones 
    
    Processing Zone: 1
    Processing Zone: 0
    
    Processing first/last row
    
    Processing Zone: 84
    ...
    
    Processing all remaining zones
    Processing Zone: 132
    Processing Zone: 4094
     merging z03331n00085c00138 
     merging z03331n00085c00138 
     merging z01968n00170c00037 
     merging z01968n00170c00037 
    ...

After this command completes, there will not be any new files. Instead, the script
re-arranges (and rewrites) any `-container.fredbin` files which have overlapping
entries. Note, because of the internal tree representation, zones 0 and 1 may be
checked multiple times.

Lastly we convert the stellar measurements in the containers (e.g.
`*-container.fredbin` files) into summarized data entries. This is accomplished
using the `rect-to-dat.py` script. There is a specific script for APASS and
SRO data. Here is some representative output data for the SRO data set above:

    ~/workspace/apass/python$ python sro_rect_to_dat.py -j8 /home/data/sro-test/*-container.fredbin 
    Processing zone z00650
    Processing zone z00688
    Processing zone z00652
    Processing zone z00705
    ...
    Processing zone z04089
    Time elapsed: 746s

The `rect-to-dat` scripts parse each zone, filters the data, averages the resulting
data, and writes the results along with some metadata to a `.dat` file. The
SRO-specific script will also generate a NetworkX Graph file (`.p`) which 
contains information on containers filled with data from multiple fields.
The output directory should look similar to the following:

    /home/data/sro-test$ ls | head -n 15
    fred-to-zone-modified-files.txt
    global.json
    z00650-border-rects.json
    z00650-container.fredbin
    z00650-contrib.txt
    z00650.dat
    z00650.fredbin
    z00650.p
    ...

Once you are satisfied with the data reduction, simply concatinate the `.dat` files.
On Linux- or Unix-based machines, you can do this:

    cat /home/data/apass/test/*.dat > output.dat
    
This final data file can be supplied to later stages of the data reduction pipeline.

## SRO Data

In the case of SRO data, the pipeline follows the same sequence as the APASS data, 
but has an additional step:

    make-zones.py
    fred-to-zone.py [-j N] *.fred
    zone-to-rect.py [-j N] *.fredbin
    fix-zone-overlaps.py [-j N]
    sro-rect-to-data.py [-j N] *-container.fredbin
    sro-merge-dat.py z*.dat
