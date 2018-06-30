#!/usr/bin/python
import argparse
import numpy as np
import numpy.lib.recfunctions as nprf
from functools import partial
import multiprocessing as mp
import os

import fred

def fredbin_add_use_data(filename):
    """Reads in a FREDBIN file, applies modifications, saves the file.
    Returns the filename upon successful modification."""

    print("Upgrading %s" % (filename))

    # Read in the FREDBIN without the
    col_names = fred.fredbin_col_names[:-1]
    col_types = fred.fredbin_col_types[:-1]
    col_fmt   = fred.fredbin_col_fmt[:-1]

    dtype={'names': col_names,'formats': col_types}
    data = np.fromfile(filename, dtype)

    append_col_names = ['use_data']
    append_col_types = ['bool']
    append_col_fmts  = ['%1i']

    num_rows = len(data)
    tmp=np.ones(num_rows)
    data = nprf.append_fields(data, append_col_names, [tmp], dtypes=append_col_types)

    fred.write_fredbin(filename, data)

    return filename

def main():

    parser = argparse.ArgumentParser(description='Parses .fred files into zone .fredbin files')
    parser.add_argument('input', nargs='+', help="Input files which will be split into zonefiles")
    parser.add_argument('-j','--jobs', type=int, help="Parallel jobs", default=4)
    parser.add_argument('--debug', default=False, action='store_true',
                        help="Run in debug mode")
    parser.set_defaults(jobs=1)

    # parse the command line arguments and start timing the script
    args = parser.parse_args()

    # load globals
    save_dir = os.path.dirname(os.path.abspath(args.input[0]))
    error_filename = save_dir + "/modify_fredbin.errorlog"

    # truncate the error log file
    with open(error_filename, 'w') as error_file:
        error_file.truncate()


    # Construct a partial to serve as the function to call in serial or
    # parallel mode below.
    fred_func = partial(fredbin_add_use_data)

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

    # write out a file containing information on the zones modified.
    mod_file = save_dir + "/" + "modified_files.txt"
    with open(mod_file, 'w') as outfile:
        for filename in results:
            outfile.write(filename + "\n")

    print("A list of modified files has been written to %s" % (mod_file))


if __name__ == "__main__":
    main()
