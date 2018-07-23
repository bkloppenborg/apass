#!/usr/bin/python

# system includes
import argparse
import glob
import re
import os
from functools import partial
import multiprocessing as mp
import numpy as np
import numpy.lib.recfunctions as nprf

# project includes
import apass
import fred

def read_contrib_file(filename):
    """Reads in a contrib file. Returns it as a numpy structured array
    with columns 'filename' and 'night_name'"""

    # read in the file as-is.
    data = []
    with open(filename, 'r') as infile:
        for line in infile:
            line = line.strip()
            night_name = find_night_name(line)
            data.append((line, night_name))

    num_rows = len(data)

    # convert data into a structured numpy array
    dtype = {'names': ['filename', 'night_name'], 'formats': ['S255', 'S7']}
    data = np.array(data, dtype=dtype)

    return data

def write_contrib_file(filename, data):
    """Writes a contrib file from a numpy structured array containing the
    fields 'filename', and 'night_name'"""

    # only the filename is written to disk
    data = data['filename']

    np.savetxt(filename, data, fmt=['%s'])

def purge_nights(night_names, save_dir, filename):
    """Loads a fredbin file, purges all data corresponding to night_name"""

    results = []

    # load the zone contribution file
    zone_id = apass.zone_from_name(filename)
    contrib_filename = save_dir + '/' + apass.name_zone_contrib_file(zone_id)
    contrib_data = read_contrib_file(contrib_filename)

    # Find common nights between the contrib file and the nights to
    # be purged. If none exist, return.
    common_nights = np.intersect1d(contrib_data['night_name'], night_names)
    if len(common_nights) == 0:
        print("No purge-able nights for zone %i" % (zone_id))
        return results

    # remove the data from the fredbin and containerized fredbin file
    fredbin_file   = save_dir + '/' + apass.name_zone_file(zone_id)
    r = purge_nights_from_file(common_nights, fredbin_file)
    results.append(r)
    container_file = save_dir + '/' + apass.name_zone_container_file(zone_id)
    r = purge_nights_from_file(common_nights, container_file)
    results.append(r)

    # now purge the nights from the contrib file
    indices = np.nonzero(np.in1d(contrib_data['night_name'], common_nights))
    contrib_data = np.delete(contrib_data, indices)
    write_contrib_file(contrib_filename, contrib_data)

def purge_nights_from_file(night_names, filename):

    results = []

    if(not os.path.isfile(filename)):
        return

    print("Processing %s" % (filename))

    # At least one common night exists. Load the data.
    data = fred.read_fredbin(filename)
    if data is None:
        return results

    # first determine if this night has been stored in this data file
    data = np.sort(data, order='night_name')

    for night_name in night_names:

        # find the upper and lower bounds for the night
        l = np.searchsorted(data['night_name'], night_name, side='left')
        r = np.searchsorted(data['night_name'], night_name, side='right')

        # delete the entries
        data = np.delete(data, slice(l, r), axis=0)

    # write the remaining data to disk
    with open(filename, 'w') as outfile:
        fred.write_fredbin(outfile, data)

    return results

def find_night_name(some_string):
    pattern = "[ns][0-9]{6}"

    night_name = None

    match = re.search(pattern, some_string)
    if match:
        l = match.start()
        r = match.end()
        night_name = some_string[l:r]

    return night_name

def main():

    parser = argparse.ArgumentParser(description='Purges data for a given night from all fredbin files')
    parser.add_argument('save_dir', help="Directory to save the output files.")
    parser.add_argument('input', nargs='+', help="Files or night names to purge")
    parser.add_argument('-j', '--jobs', type=int, help="Parallel jobs", default=4)
    parser.add_argument('--debug', default=False, action='store_true',
                        help="Run in debug mode")

    parser.set_defaults(jobs=1)
    args = parser.parse_args()

    # Determine the night names that should be deleted
    night_names = []
    for item in args.input:
        night_name = find_night_name(item)

        if night_name != None:
            night_names.append(night_name)
        else:
            print("Skipping %s, could not determine night name" % (item))

    # Find all of the fredbin files.
    fredbin_files = glob.glob(args.save_dir + "*-raw.fredbin")

    # construct the purge function and run the deletion process
    purge_func = partial(purge_nights, night_names, args.save_dir)
    results = []
    if args.debug:
        for filename in fredbin_files:
            r = purge_func(filename)
            results.append(r)
    else:
        pool = mp.Pool(args.jobs)

        # farm out the work and wait for the result
        results = pool.imap(purge_func, fredbin_files)
        pool.close()
        pool.join()

if __name__ == "__main__":
    main()
