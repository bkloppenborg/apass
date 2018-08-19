#!/usr/bin/python

import argparse
import numpy as np
import multiprocessing as mp
from functools import partial

# local includes
import fred
import apass

from dat import apass_phot_names as apass_phot_names # named values
from dat import apass_filter_ids as apass_filter_ids # numerical values

def summary_dict_cols():
    """Supplies the column names, column formats, and column types
    for the summary dictionary."""

    col_names = ['filename', 'field_id']
    col_names.extend(apass_phot_names)
    col_fmts = ['%50s', '%25s']
    col_fmts.extend(['%5i'] * len(apass_phot_names))
    col_types = ['S50', 'S25']
    col_types.extend(['int'] * len(apass_phot_names))

    return col_names, col_fmts, col_types

def init_summary_dict(filename, fred_type="apass"):
    """Creates and initializes a summary dictionary"""

    output = dict()
    output['filename'] = filename
    output['ra']       = 0
    output['dec']      = 0
    output['field_id'] = ""

    for filter_name in apass_phot_names:
        output[filter_name] = 0

    return output

def output_dict_cols():
    """Supplies the colmn names, column formats, and column types
    for the output dictionary."""

    col_names = ['field_id', 'ra', 'dec']
    col_names.extend(apass_phot_names)
    col_fmts = ['%25s', '%11.6f', '%11.6f']
    col_fmts.extend(['%5i'] * len(apass_phot_names))
    col_types = ['S25', 'float64', 'float64']
    col_types.extend(['int'] * len(apass_phot_names))

    return col_names, col_fmts, col_types

def init_output_dict():
    """Creates and initializes an output dictionary"""

    output = dict()
    output['ra'] = 0
    output['dec'] = 0
    output['field_id'] = 0

    for filter_name in apass_phot_names:
        output[filter_name] = 0

    return output

def summarize_fred(filename):
    """Reads in a FRED/DTFRED file, generates some basic statistics on a per-filter basis"""

    print("Processing %s" % (filename))

    output = []

    data = None
    try:
        data = fred.read_fred(filename)
    except:
        pass

    # bipass empty files
    if data is None:
        return output

    # find unique observations of fields in each filter
    field_ids = set(data['field_id'])

    # for each field, summarize the number of observations in each filter
    for field_id in field_ids:
        # init a summary structure, start populating the data
        t_summary = init_summary_dict(filename)
        t_summary['field_id'] = field_id

        # extract data for this field
        indices = np.where(data['field_id'] == field_id)
        field_data = data[indices]

        # for each filter, determine the number of unique observations
        for filter_id in apass_filter_ids:

            # extract the data for this filter
            indices = np.where(field_data['filter_id'] == filter_id)
            filter_data = field_data[indices]

            # Get the unique JDs for this filter. This is a proxy for the
            # number of distinct observations.
            unique_sets = set(filter_data['set'])

            # save the number of observations
            filter_name = apass.filter_name_from_id(filter_id)
            t_summary[filter_name] = len(unique_sets)

        output.append(t_summary)

    return output

def dicts_to_ndarray(dicts, col_names, col_types):
    """Converts a dictionary to a structured numpy array in a .dat-friendly format"""

    # convert each dict into a list in the correct order
    output = []
    for d in dicts:
        # convert the dictionary to a list in the right order
        t = []
        for k in col_names:
            t.append(d[k])

        # grow the output
        output.append(tuple(t))

    # convert the output to a numpy named array
    dtype={'names': col_names, 'formats': col_types}
    output = np.asarray(output, dtype=dtype)
    return output

def lookup_field_centers(field_id, field_centers):
    """Finds the matching field center in a SORTED numpy structured array
    containing the field center data."""

    num_fields = len(field_centers)

    # find the matching field. This should be unique.
    idx = np.searchsorted(field_centers['field_id'], field_id)

    ra  = 0.0
    dec = 0.0

    if idx < num_fields and field_centers['field_id'][idx] == field_id:
        ra  = field_centers[idx]['ra']
        dec = field_centers[idx]['dec']

    return ra, dec

def save_night_summary(filename, data):
    """Saves information about each night to a file."""
    pass


def save_field_summary(filename, dicts, field_centers):
    """Saves information on a given field to a file"""
    output = []

    # convert the summary dictionaries into a numpy structured array
    col_names, col_fmts, col_types = summary_dict_cols()
    data = dicts_to_ndarray(dicts, col_names, col_types)

    # find all of the unique field entries, add them up
    field_ids = set(data['field_id'])
    for field_id in field_ids:

        ra, dec = lookup_field_centers(field_id, field_centers)
        field_output             = init_output_dict()
        field_output['field_id'] = field_id
        field_output['ra']       = ra
        field_output['dec']      = dec

        indices = np.where(data['field_id'] == field_id)
        t_data = data[indices]

        for key in apass_phot_names:
            field_output[key] = np.sum(t_data[key])

        output.append(field_output)

    # convert the output dictionary to an ndarray
    col_names, col_fmts, col_types = output_dict_cols()
    data = dicts_to_ndarray(output, col_names, col_types)

    # construct a header
    header =  "APASS FRED summary file output. Format is as follows: \n"
    header += "Column names: " + ','.join(col_names) + "\n"
    header += "Column types: " + ','.join(col_types) + "\n"
    header += "Column formats: " + ','.join(col_fmts) + "\n"

    # write the results to a file.
    np.savetxt(filename, data, fmt=col_fmts, header=header)

def read_field_centers(filename):

    col_names = ['field_id', 'ra', 'dec']
    col_types = ['S25', 'float64', 'float64']

    dtype={'names': col_names,'formats': col_types}
    data = np.loadtxt(filename, dtype=dtype)

    data = np.sort(data, order=['field_id'])

    return data

def main():

    parser = argparse.ArgumentParser(description='Prints statistical information on FRED/DTFRED files')
    parser.add_argument('field_center_file', type=str, help="File containing field center descriptions")
    parser.add_argument('input', nargs='+', help="Input files to summarize")
    parser.add_argument('-j','--jobs', type=int, help="Parallel jobs", default=4)
    parser.add_argument('--debug', default=False, action='store_true',
                        help="Run in debug mode")
    parser.set_defaults(jobs=1)

    # parse the command line arguments and start timing the script
    args = parser.parse_args()
    field_centers = read_field_centers(args.field_center_file)

    run_func = partial(summarize_fred)

    results = []
    if args.debug:
        for filename in args.input:
            r = run_func(filename)
            results.extend(r)
    else:
        pool = mp.Pool(args.jobs)

        # farm out the jobs and wait for the result
        pool_result = pool.imap(run_func, args.input)
        pool.close()
        pool.join()

        for r in pool_result:
            results.extend(r)


    save_night_summary('fred_night_summary.txt', results)
    save_field_summary('fred_field_summary.txt', results, field_centers)

if __name__ == "__main__":
    main()
