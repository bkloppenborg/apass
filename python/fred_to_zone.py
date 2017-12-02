#!/bin/python
import argparse
import numpy as np
from numpy import cos
from math import pi
import inspect
import glob

# parallel processing
import multiprocessing as mp

# APASS-specific things
from quadtree import *
from quadtree_types import *
import apass
from apass import *

import sys, os
sys.path.append(os.path.join(sys.path[0],'modules', 'FileLock', 'filelock'))
from filelock import FileLock

def build_data_dict(filename):
    """Creates a dictionary which maps the data in the specified file to specific
    zone IDs."""

    data_dict = {}

    global tree_file

    # restore the tree. make-zones.py writes out leaves of type IDLeaf
    tree = QuadTreeNode.from_file(tree_file, leafClass=IDLeaf)

    # read the fred file into a numpy array
    print("Processing FRED file " + filename)
    try:
        data = read_fred(filename)
    except ValueError:
        print("ERROR: File %s has a bad value and was not parsed" % (filename))
        return None
    except IndexError:
        print("ERROR: File %s is missing data and was not parsed" % (filename))
        return None
    except:
        print("ERROR: File %s has an unknown error and was not parsed" % (filename))
        return None

    # process every file, inserting it into a .dat file. Keep track of any
    # files that were opened in the open_files set
    open_files = set()
    for datum in np.nditer(data):
        # we are about to modify the datum, make a copy
        datum = datum.copy()

        # pull out the RA and DEC
        [ra, dec] = get_coords(datum)
        zone_id = tree.insert(ra, dec, None)

        if zone_id not in data_dict.keys():
            data_dict[zone_id] = list()

        # update the zone ID and append it to the dictionary
        datum['zone_id'] = zone_id
        data_dict[zone_id].append(datum)

    return data_dict


def add_fred(filename):
    """Processes an APASS FRED file into zones using data contained in the tree"""
    data_dict = build_data_dict(filename)
    impacted_zones = []

    if data_dict is not None:
        # Write out the data being sure to lock all related files prior to opening
        for zone_id, data in data_dict.iteritems():
            zone_filename = apass_save_dir + '/' + name_zone_file(zone_id)
            contrib_filename = apass_save_dir + '/' + name_zone_contrib_file(zone_id)

            with FileLock(zone_filename, timeout=100, delay=0.05):
                with open(zone_filename, 'a+b') as outfile:
                    for datum in data:
                        outfile.write(datum)

                with open(contrib_filename, 'a+') as outfile:
                    outfile.write(filename + "\n")

            impacted_zones.append(zone_id)

    return impacted_zones

def remove_fred(filename):
    """Removes the data found in the specified file."""

    data_dict = build_data_dict(filename)
    impacted_zones = []

    if data_dict is not None:
        # Write out the data being sure to lock all related files prior to opening
        for zone_id, data in data_dict.iteritems():

            # print out a message to inform the user of progress.
            print("Checking zone %i" % (zone_id))
            data_removed = False

            zone_filename = apass_save_dir + '/' + name_zone_file(zone_id)
            contrib_filename = apass_save_dir + '/' + name_zone_contrib_file(zone_id)

            with FileLock(zone_filename):
                zone_data = apass.read_fredbin(zone_filename)

                zone_data = np.sort(zone_data, order=['ra', 'dec'])
                num_data = len(zone_data)

                # find the indices of the corresponding data points
                # ensure that the data match *exactly)
                indices = []
                for datum in data:
                    index = np.searchsorted(zone_data, datum)
                    if index >= num_data:
                        continue

                    same_point = apass.compare_fred_data(datum, zone_data[index])
                    if same_point:
                        indices.append(index)
                        data_removed = True

                print("Deleting %i entries from zone %i" % (len(indices), zone_id))
                zone_data = np.delete(zone_data, indices)

                # write the file
                with open(zone_filename, 'w+b') as outfile:
                    for datum in zone_data:
                        outfile.write(datum)

            if data_removed:
                impacted_zones.append(zone_id)


    return impacted_zones

def main():

    parser = argparse.ArgumentParser(description='Parses .fred files into zone .fredbin files')
    #parser.add_argument('outdir', help="Directory into which .fredbin files will be generated")
    parser.add_argument('input', nargs='+', help="Input files which will be split into zonefiles")
    parser.add_argument('-j','--jobs', type=int, help="Parallel jobs", default=4)
    parser.add_argument('--debug', default=False, action='store_true',
                        help="Run in debug mode")
    parser.add_argument('--remove', default=False, action='store_true')
    parser.set_defaults(jobs=1)

    args = parser.parse_args()
    outdir = apass_save_dir

    # load the global tree
    global tree_file
    tree_file = apass_save_dir + "/global.json"

    # override the function we use if we are removing data
    fred_func = add_fred
    if args.remove:
        fred_func = remove_fred

    # set up the pool and launch the function
    results = []
    if args.debug:
        for filename in args.input:
            r = fred_func(filename)
            results.extend(r)
    else:
        pool = mp.Pool(args.jobs)

        # farm out the jobs and wait for the result
        pool_result = pool.imap(fred_func, args.input)
        pool.close()
        pool.join()

        for r in pool_result:
            results.extend(r)

    # determine the zone files impacted by this operation:
    zone_names = []
    zone_ids = sorted(set(results))
    for zone_id in zone_ids:
        filename = apass.name_zone_file(zone_id)
        zone_names.append(filename)

    # write out a file containing information on the zones modified.
    mod_file = apass_save_dir + "fred-to-zone-modified-files.txt"
    with open(mod_file, 'w') as outfile:
        for filename in zone_names:
            outfile.write(apass_save_dir + filename + "\n")

    print("A list of modified files has been written to %s" % (mod_file))


if __name__ == "__main__":
    main()
