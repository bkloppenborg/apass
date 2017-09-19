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

These packages can be installed on a Ubuntu system with the following commands:

    sudo apt install python-pip python-numpy python-matplotlib
    pip install numpy argparse matplotlib python-tk
    pip install nxgraph pydotplus


## Sample reduction run

Below is a worked example of running the data reduction pipeline on a subset of
the APASS DR9 FRED files. Sample output from each programs is listed in addition
to contents found in the output directories.

For the purposes of this sample, the following directories are used:

    APASS Code:           ~/workspace/apass/python 
    FRED directory:       /home/data/test-freds
    APASS save directory: /home/data/apass-test
    
Create zones using `make-zones.py`

    ~/workspace/apass/python$ python make-zones.py 
    Created 3904 zones
    
This will result in a `global.json` file in the APASS save directory:

    ~/workspace/apass/python$ ls /home/data/apass-test/
    global.json
   
Next we parse each FRED file and split the data into different zones. To run
this process, we use the `fred-to-zone.py` script with `.fred` files as input:

    ~/workspace/apass/python$ python fred-to-zone.py /home/data/test-freds/*.fred -j 8
    Processing FRED file /home/data/test-freds/130112.fred
    Processing FRED file /home/data/test-freds/130105.fred
    Processing FRED file /home/data/test-freds/130103.fred
    Processing FRED file /home/data/test-freds/130107.fred
    Processing FRED file /home/data/test-freds/130114.fred
    Processing FRED file /home/data/test-freds/130127.fred
    Processing FRED file /home/data/test-freds/130202.fred
    Processing FRED file /home/data/test-freds/130125.fred
    ...
    
This command will generate a series of `.fredbin` files along with a
`*-contrib.txt` file as seen in the output directory below:

    ~/workspace/apass/python$ ls /home/data/apass-test/ | head -n 5
    global.json
    z00001-contrib.txt
    z00001.fredbin
    z00008-contrib.txt
    z00008.fredbin
  
  
The `.fredbin` files are a binary representation of the `.fred` files with a few
additonal columns (for bookkeeping zone, node, and container IDs) appended to
the end. If you wish to see the content of the `.fredbin`, there is a `dump-fredbin.py`
script included with the data reduction pipeline.

Next we parse the zone data and segment the measurements into separate bins, one
for each star. Internally, the bins are represented as bounding rectangles whose
bounds encompass all of the measurements for a star, plus a small (few arcsecond)
border to serve as padding. This stage of the pipeline is accomplished using the
`zone-to-rects.py` script:

    ~/workspace/apass/python$ python zone-to-rects.py /home/data/apass-test/*.fredbin -j 8
    Processing '/home/data/apass-test/z00001.fredbin' which has 231 data points 
    Processing '/home/data/apass-test/z00056.fredbin' which has 3342 data points 
    Processing '/home/data/apass-test/z00130.fredbin' which has 2714 data points 
    Processing '/home/data/apass-test/z00171.fredbin' which has 634 data points 
    Processing '/home/data/apass-test/z00212.fredbin' which has 6227 data points 
    Processing '/home/data/apass-test/z00255.fredbin' which has 20577 data points 
    Processing '/home/data/apass-test/z00418.fredbin' which has 27033 data points
    ...
    
After this command executes, you will find directory contents similar to the following:
    
    ~/workspace/apass/python$ ls /home/data/apass-test/ | head -n 10
    global.json
    z00001-border-rects.json
    z00001-container.fredbin
    z00001-contrib.txt
    z00001.fredbin
    z00001-zone.json
    z00008-border-rects.json
    z00008-container.fredbin
    z00008-contrib.txt
    z00008.fredbin

The `*-container.fredbin` files contain the same data as the previous `.fredbin`
files, but are now sorted/grouped by star and contain non-zero entries for the
node and container IDs. The script also generates a series of
`*-border-rects.json` files which describe any rectangles that might span more
than one zone.
    
Next we check for any overlaps between zones using `fix-zone-overlaps.py`. At
present this step does not have a parallel execution option. It automatically
parses data found in the APASS save directory, so simply run the script as
follows:

    ~/workspace/apass/python fix-zone-overlaps.py -j 8 
    Processing Zone: 1
     checking for overlaps with z00001n00171c00060
      merging z00001n00171c00060 
      updating border rect file for zone 1
     checking for overlaps with z00001n00170c00062
     checking for overlaps with z00001n00170c00063
    Processing Zone: 1
     checking for overlaps with z00001n00170c00062
     checking for overlaps with z00001n00170c00063
    Processing Zone: 8
    ...

After this command completes, there will not be any new files. Instead, the script
re-arranges (and rewrites) any `-container.fredbin` files which have overlapping
entries. Note, because of the internal tree representation, zones 0 and 1 may be
checked multiple times.

Lastly we convert the stellar measurements in the containers (e.g.
`*-container.fredbin` files) into summarized data entries. This is accomplished
using the `rect-to-data.py` script as follows:

    ~/workspace/apass/python$ python apass-rect-to-data.py /home/data/apass-test/*-container.fredbin -j 8
    Processing zone z00055
    Processing zone z00001
    Processing zone z00129
    Processing zone z00169
    Processing zone z00210
    Processing zone z00252
    Processing zone z00415
    Processing zone z00481
    ...

It will generate a series of `.dat` files that contain averaged RA, DEC, one
magnitude for each APASS photometric filter, and all corresponding
uncertainties:
    
    ~/workspace/apass/python$ ls /home/data/apass-test/ | head -n 10
    global.json
    z00001-border-rects.json
    z00001-container.fredbin
    z00001-contrib.txt
    z00001.dat
    z00001.fredbin
    z00001-zone.json
    z00008-border-rects.json
    z00008-container.fredbin
    z00008-contrib.txt

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
    sro-rect-to-data.py [-j N] *-container.fredbinA
    sro-merge-dat.py
