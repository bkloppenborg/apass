#!/usr/bin/python

import argparse
import numpy as np
import multiprocessing as mp
from functools import partial

# local includes
import fred

def init_summary_dict(filename):

    output = dict()
    output['filename'] = filename
    output['entries']  = 0
    output['nights']   = dict()

    return output

def summarize_fred(filename):
    """Reads in a fredbin file, prints out some basic statistics"""

    print("Processing %s" % (filename))

    output = init_summary_dict(filename)

    data = None
    try:
        data = fred.read_fredbin(filename)
    except:
        pass

    # bipass empty files
    if data is None:
        return output

    # populate the dictionary with useful information
    data = np.sort(data, order='night_name')
    output['entries'] = len(data)

    # find the unique named nights in the data
    night_names = list(set(data['night_name']))
    for night_name in night_names:

        # find the upper and lower bounds for the night
        l = np.searchsorted(data['night_name'], night_name, side='left')
        r = np.searchsorted(data['night_name'], night_name, side='right')

        output['nights'][night_name] = r - l

    return output

def save_night_summary(filename, data):

    # extract the night information from the results
    night_data = dict()

    for d in data:
        temp = d['nights']
        keys = temp.keys()

        for k in keys:
            if k in night_data.keys():
                night_data[k] += temp[k]
            else:
                night_data[k]  = temp[k]

    # write the results to a file
    with open(filename, 'w') as outfile:
        keys = night_data.keys()

        for k in keys:
            k_value = str(night_data[k])
            outfile.write(k + " " + k_value + "\n")


def save_zone_summary(filename, data):

    with open(filename, 'w') as outfile:
        for d in data:
            filename = d['filename']
            entries  = d['entries']

            outfile.write(filename + " " + str(entries) + "\n")


def main():

    parser = argparse.ArgumentParser(description='Prints statistical information on FREDBIN files')
    parser.add_argument('input', nargs='+', help="Input files which will be split into zonefiles")
    parser.add_argument('-j','--jobs', type=int, help="Parallel jobs", default=4)
    parser.add_argument('--debug', default=False, action='store_true',
                        help="Run in debug mode")
    parser.set_defaults(jobs=1)

    # parse the command line arguments and start timing the script
    args = parser.parse_args()

    run_func = partial(summarize_fred)

    results = []
    if args.debug:
        for filename in args.input:
            r = run_func(filename)
            results.append(r)
    else:
        pool = mp.Pool(args.jobs)

        # farm out the jobs and wait for the result
        pool_result = pool.imap(run_func, args.input)
        pool.close()
        pool.join()

        for r in pool_result:
            results.append(r)

    save_night_summary('night_summary.txt', results)
    save_zone_summary('zone_summary.txt', results)

if __name__ == "__main__":
    main()
