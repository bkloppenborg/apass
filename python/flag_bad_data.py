#!/usr/bin/python
import argparse
import numpy as np
import os
import traceback
import multiprocessing as mp
from functools import partial

# lockfile
import sys, os
sys.path.append(os.path.join(sys.path[0],'modules', 'FileLock', 'filelock'))
from filelock import FileLock

# local includes
import badfiles
import fred

def flag_bad_data(bad_night_file, bad_field_file, filename):
    """Sets the 'use_data' entry to FALSE for any data taken on a known bad night or
    bad night-field combination."""

    global error_filename

    try:

        print("Processing %s" % (filename))

        # read in the bad night/field files
        bad_night_names = badfiles.read_bad_nights(bad_night_file)
        bad_fields = badfiles.read_bad_night_fields(bad_field_file)

        # read in the data and sort it by night
        data = fred.read_fredbin(filename)

        # Any given filter should ONLY set 'use_data' flags to False to avoid
        # impacting other filters.
        data = np.sort(data, order='night_name')
        data = filter_bad_nights(data, bad_night_names)
        data = np.sort(data, order='night')
        data = filter_bad_night_fields(data, bad_fields)
        data = filter_non_photometric_nights(data)

        # restore the original order of the data
        data = np.sort(data, order=['zone_id', 'node_id', 'container_id'])

        # write the data
        outfile = open(filename, 'w')
        fred.write_fredbin(outfile, data)

    except:
        message = "ERROR: Failed to flag %s. Re-run in debug mode\n" % (filename)

        tb = traceback.format_exc()

        print(message)
        with FileLock(error_filename, timeout=100, delay=0.05):
            with open(error_filename, 'a') as error_file:
                error_file.write(message + "\n" + str(tb) + "\n")

def filter_non_photometric_nights(data):
    """Sets the 'use_data' flag to False for nights that are identified as
    being non-photometric in the raw FRED fields. Returns the modified data array"""

    indexes = np.where(data['flag1'] == 1)
    data['use_data'][indexes] = False

    return data


def filter_bad_nights(data, bad_night_names):
    """Sets the 'use_data' flag to False for nights that are identified as
    bad nights.
    Input:
    data - numpy array loaded using fred.read_fredbin
    bad_nights - numpy array loaded using badfiles.read_bad_nights

    Returns:
    modified numpy array in fredbin format
    """

    # iterate over the bad night entries, flagging data as required
    for i in range(0, len(bad_night_names)):
        # extract the bad (numeric) night and corresponding night name
        night_name  = bad_night_names['night_name'][i]  # string

        # find upper and lower bounds for the night
        l = np.searchsorted(data['night_name'], night_name, side='left')
        u = np.searchsorted(data['night_name'], night_name, side='right')

        # iterate over each candidate and flag if the night name
        for j in range(l, u):
            data['use_data'][j] = False

    return data

def filter_bad_night_fields(data, bad_nights_fields):
    """Sets the 'use_data' flag to False for fields on specific nights that
    have been identified as bad.  Returns the modified data array"""

    for i in range(0, len(bad_nights_fields)):
        # extract the bad (numeric) night and corresponding field name
        night = bad_nights_fields['night'][i] # integer
        field_id = bad_nights_fields['field_id'][i] # string

        # find upper and lower bounds for the night
        l = np.searchsorted(data['night'], night, side='left')
        u = np.searchsorted(data['night'], night, side='right')

        # iterate over each candidate night and flag if the field
        # name is a precise match.
        for j in range(l, u):
            if data['field_id'][j] == field_id:
                data['use_data'][j] = False

    return data

def main():

    parser = argparse.ArgumentParser(
        description="Converts a containerized zone into APASS photometric output")
    parser.add_argument('bad_night_file', help="File containing nights that should be excluded")
    parser.add_argument('bad_field_file', help="File containing night-fields that should be excluded")
    parser.add_argument('input', nargs="+",
                        help="FREDBIN files to flag.")
    parser.add_argument('-j','--jobs', type=int, help="Parallel jobs", default=4)
    parser.add_argument('--debug', default=False, action='store_true',
                        help="Run in debug mode")
    args = parser.parse_args()

    # specify globals
    global error_filename

    # Extract the save directory from the 1st file input
    save_dir = os.path.dirname(os.path.realpath(args.input[0])) + "/"

    # clear out the error file
    error_filename = save_dir + "error_flag_bad_data.txt"
    try:
        os.remove(error_filename)
    except:
        pass

    # Configure the flag bad data function
    flag_func = partial(flag_bad_data, args.bad_night_file, args.bad_field_file)

    # run in debug mode
    if args.debug:
        for filename in args.input:
            flag_func(filename)
    # run in production mode
    else:
        pool = mp.Pool(args.jobs)
        result = pool.imap(flag_func, args.input)
        pool.close()
        pool.join()




if __name__ == "__main__":
    main()
