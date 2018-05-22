#!/bin/python

# system includes
import argparse
import numpy as np
from numpy import cos
from math import pi
import inspect
import glob
from functools import partial
import time
import traceback
import datetime

# parallel processing
import multiprocessing as mp

# APASS-specific things
from quadtree import *
from quadtree_types import *
from apass import name_zone_file, name_zone_contrib_file, name_zone_file
from apass import get_coords, get_num_zones

# File I/O
from fred import read_fred, read_fredbin, compare_fred_data

import sys, os
sys.path.append(os.path.join(sys.path[0],'modules', 'FileLock', 'filelock'))
from filelock import FileLock

def make_data_dict():
    """Constructs a data dictionary consisting of the following fields:
       num_fred_data - integer, number of values read in from FRED files
       XXXXXX        - zone ID as key with data /mnt/astro/astro/home/welch/surveys/APASS/hphot_and_hepp/fixedk2_dtfred/s160201.fixedk2.dtfredto be inserted into that zone as value
    """

    out = dict()
    out['num_fred_data'] = 0
    num_zones = get_num_zones()
    for i in range(0, num_zones):
        out[i] = list()

    return out

def build_data_dict(filename):
    """Creates a dictionary which maps the data in the specified file to specific
    zone IDs."""

    global tree_file
    global error_filename

    # create the data dictionary which will store data values.
    data_dict = make_data_dict()

    flag_read_error = False

    # restore the tree. make-zones.py writes out leaves of type IDLeaf
    tree = QuadTreeNode.from_file(tree_file, leafClass=IDLeaf)

    # read the fred file into a numpy array
    try:
        data = read_fred(filename)
    except ValueError:
        print("ERROR: File %s has a bad value and was not parsed" % (filename))
        flag_read_error = True
    except IndexError:
        print("ERROR: File %s is missing data and was not parsed" % (filename))
        flag_read_error = True
    except:
        print("ERROR: File %s has an unknown error and was not parsed" % (filename))
        flag_read_error = True

    if flag_read_error:
        with FileLock(error_filename, timeout=100, delay=0.05):
            with open(error_filename, 'a') as error_file:
                error_file.write("ERROR: Cannot parse %s" % (filename))
        return None

    # update the number of data points read.
    data_dict['num_fred_data'] = len(data)

    # Map each data point to a corresponding zone file by performing
    # a mock-insert into the (local copy) of the global tree and storing
    # the data point within

    open_files = set()
    for datum in np.nditer(data):
        # we are about to modify the datum, make a copy
        datum = datum.copy()

        # pull out the RA and DEC and perform a mock insert into the global tree
        [ra, dec] = get_coords(datum)
        zone_id = tree.insert(ra, dec, None)

        # update the datum with its new zone ID
        datum['zone_id'] = zone_id

        # move the data into the dictionary
        data_dict[zone_id].append(datum)

    return data_dict

def write_mapping_info(save_dir, filename, data_dict, mode="add"):
    """
    Valid modes are "add" or "remove"
    """

    num_zones = get_num_zones()
    date = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    # sign to indicate addition or removal of data
    sign = 1
    if mode == "remove":
        sign = -1

    # write out the input to zone assignment information
    mapping_data = list()
    mapping_data.append(date)
    mapping_data.append(filename)
    mapping_data.append(data_dict['num_fred_data'])
    for zone_id in range(0, num_zones):
        num_data = 0
        zone_data = data_dict[zone_id]
        if zone_data is not None:
            num_data = len(zone_data)

        num_data *= sign
        mapping_data.append(num_data)

    # convert the data to strings
    mapping_data = map(str, mapping_data)

    # write the data to file
    zone_mapping_file = save_dir + '/' + 'zone_mapping.log'
    with FileLock(zone_mapping_file, timeout=100, delay=0.05):
        with open(zone_mapping_file, 'a+') as outfile:
            outfile.write(','.join(mapping_data) + "\n")

def add_fred(save_dir, filename):
    """Processes an APASS FRED file into zones using data contained in the tree"""
    print("Processing FRED file " + filename)

    impacted_zones = []
    data_dict = build_data_dict(filename)

    # if there isn't any data, bail out early.
    if data_dict is None:
        return impacted_zones

    # write the input data to zone mapping information to a file
    write_mapping_info(save_dir, filename, data_dict, mode="add")

    # remove any data that is not for a zone
    del data_dict['num_fred_data']

    # Write out the data being sure to lock all related files prior to opening
    for zone_id, data in data_dict.iteritems():

        # skip zones with no data
        if len(data) == 0:
            continue

        zone_filename    = save_dir + '/' + name_zone_file(zone_id)
        contrib_filename = save_dir + '/' + name_zone_contrib_file(zone_id)

        with FileLock(zone_filename, timeout=100, delay=0.05):
            with open(zone_filename, 'a+b') as outfile:
                for datum in data:
                    outfile.write(datum)

            with open(contrib_filename, 'a+') as outfile:
                outfile.write(filename + "\n")

        impacted_zones.append(zone_id)

    print("Completed FRED file " + filename)

    return impacted_zones

def remove_fred(save_dir, filename):
    """Removes the data found in the specified file."""

    print("Processing FRED file " + filename)
    impacted_zones = []
    data_dict = build_data_dict(filename)

    # if there isn't any data, bail out early.
    if data_dict is None:
        return impacted_zones

    # write the input data to zone mapping information to a file
    write_mapping_info(save_dir, filename, data_dict, mode="remove")

    # remove any data that is not for a zone
    del data_dict['num_fred_data']

    if data_dict is not None:
        # Write out the data being sure to lock all related files prior to opening
        for zone_id, data in data_dict.iteritems():

            if len(data) == 0:
                continue

            # print out a message to inform the user of progress.
            #print("Checking zone %i" % (zone_id))
            data_removed = False

            zone_filename    = save_dir + '/' + name_zone_file(zone_id)
            contrib_filename = save_dir + '/' + name_zone_contrib_file(zone_id)

            with FileLock(zone_filename):
                try:
                    zone_data = read_fredbin(zone_filename)
                except:
                    print("Warning, could not find zone %i" % (zone_id))
                    continue

                zone_data = np.sort(zone_data, order=['ra', 'dec'])
                num_data = len(zone_data)

                # find the indices of the corresponding data points
                # ensure that the data match *exactly)
                indices = []
                for datum in data:
                    index = np.searchsorted(zone_data, datum)
                    if index >= num_data:
                        continue

                    same_point = compare_fred_data(datum, zone_data[index])
                    if same_point:
                        indices.append(index)
                        data_removed = True

                #print("Deleting %i entries from zone %i" % (len(indices), zone_id))
                zone_data = np.delete(zone_data, indices)

                # write the file
                with open(zone_filename, 'w+b') as outfile:
                    for datum in zone_data:
                        outfile.write(datum)

            if data_removed:
                impacted_zones.append(zone_id)

    print("Completed FRED file " + filename)

    return impacted_zones

def main():

    parser = argparse.ArgumentParser(description='Parses .fred files into zone .fredbin files')
    parser.add_argument('input', nargs='+', help="Input files which will be split into zonefiles")
    parser.add_argument('save_dir', help="Directory to save the output files.")
    parser.add_argument('-j','--jobs', type=int, help="Parallel jobs", default=4)
    parser.add_argument('--debug', default=False, action='store_true',
                        help="Run in debug mode")
    parser.add_argument('--remove', default=False, action='store_true')
    parser.set_defaults(jobs=1)

    # parse the command line arguments and start timing the script
    args = parser.parse_args()
    start = time.time()

    # load the global tree
    global tree_file
    tree_file = args.save_dir + "/global.json"

    global error_filename
    error_filename = args.save_dir + "/fred_to_zone.errorlog"

    # truncate the error log file
    with open(error_filename, 'w') as error_file:
        error_file.truncate()

    # ensure the zone mapping file exists and is populated with a header
    zone_mapping_file = args.save_dir + '/' + 'zone_mapping.log'
    if not os.path.isfile(zone_mapping_file):
        num_zones = get_num_zones()
        header = list()
        header.append('date')
        header.append('fred_filename')
        header.append('fred_num_values')
        for i in range(0, num_zones):
            header.append(i)

        header = map(str, header)
        with open(zone_mapping_file, 'w') as outfile:
            outfile.write(','.join(header) + "\n")

    # Construct a partial to serve as the function to call in serial or
    # parallel mode below.
    fred_func = partial(add_fred, args.save_dir)
    if args.remove:
        fred_func = partial(remove_fred, args.save_dir)

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
        filename = name_zone_file(zone_id)
        zone_names.append(filename)

    # write out a file containing information on the zones modified.
    mod_file = args.save_dir + "fred-to-zone-modified-files.txt"
    with open(mod_file, 'w') as outfile:
        for filename in zone_names:
            outfile.write(args.save_dir + '/' + filename + "\n")

    print("A list of modified files has been written to %s" % (mod_file))

    end = time.time()
    print("Time elapsed: %is" % (int(end - start)))

if __name__ == "__main__":
    main()
