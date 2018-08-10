#!/usr/bin/python

import argparse
import numpy
import multiprocessing as mp
from functools import partial

# local includes
import fred

def init_summary_dict(filename):

    output = dict()
    output['filename'] = filename
    output['entries'] = 0

    return output

def print_summary(s):

    filename = s['filename']
    entries  = s['entries']

    print("%s %i" % (filename, entries))

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
    output['entries'] = len(data)

    return output

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

    # print the result to console
    for r in results:
        print_summary(r)

if __name__ == "__main__":
    main()
