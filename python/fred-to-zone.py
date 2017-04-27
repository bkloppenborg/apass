#!/bin/python
import argparse
import numpy as np
from numpy import cos
from math import pi
import inspect
import glob
from enum import Enum

# parallel processing
import multiprocessing as mp

# APASS-specific things
from quadtree import *
from quadtree_types import *
from apass import *

import sys, os
sys.path.append(os.path.join(sys.path[0],'modules', 'FileLock', 'filelock'))
from filelock import FileLock

def fred_to_zone_func(filename):
    """Processes an APASS FRED file into zones using data contained in the tree"""

    zone_dict = {}

    global tree_file

    # restore the tree. make-zones.py writes out leaves of type IDLeaf
    tree = QuadTreeNode.from_file(tree_file, leafClass=IDLeaf)

    # read the fred file into a numpy array
    print("Processing FRED file " + filename)
    data = read_fred(filename)

    # process every file, inserting it into a .dat file. Keep track of any
    # files that were opened in the open_files set
    open_files = set()
    for datum in np.nditer(data):
        # we are about to modify the datum, make a copy
        datum = datum.copy()

        # pull out the RA and DEC
        [ra, dec] = get_coords(datum)
        zone_id = tree.insert(ra, dec, None)

        if zone_id not in zone_dict.keys():
            zone_dict[zone_id] = list()

        # update the zone ID and append it to the dictionary
        datum['zone_id'] = zone_id
        zone_dict[zone_id].append(datum)


    # Write out the data being sure to lock all related files prior to opening
    for zone_id, data in zone_dict.iteritems():
        zone_filename = apass_save_dir + '/' + name_zone_file(zone_id)
        contrib_filename = apass_save_dir + '/' + name_zone_contrib_file(zone_id)

        with FileLock(zone_filename):
            with open(zone_filename, 'a+b') as outfile:
                for datum in data:
                    outfile.write(datum)

            with open(contrib_filename, 'a+') as outfile:
                outfile.write(filename + "\n")

def main():

    parser = argparse.ArgumentParser(description='Parses .fred files into zone .fredbin files')
    #parser.add_argument('outdir', help="Directory into which .fredbin files will be generated")
    parser.add_argument('input', nargs='+', help="Input files which will be split into zonefiles")
    parser.add_argument('-j','--jobs', type=int, help="Parallel jobs", default=4)
    parser.add_argument('--debug', default=False, action='store_true',
                        help="Run in debug mode")
    parser.set_defaults(jobs=1)

    args = parser.parse_args()
    outdir = apass_save_dir

    # load the global tree
    global tree_file
    tree_file = apass_save_dir + "/global.json"

    # set up the pool
    if args.debug:
        for filename in args.input:
            fred_to_zone_func(filename)

    else:
        pool = mp.Pool(args.jobs)

        # farm out the jobs and wait for the result
        result = pool.map_async(fred_to_zone_func, args.input)
        result.get()

        # wait for the pool to complete
        pool.close()
        pool.join()


if __name__ == "__main__":
    main()
