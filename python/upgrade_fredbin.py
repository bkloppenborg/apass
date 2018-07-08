#!/usr/bin/python
import argparse
import numpy as np
import multiprocessing as mp

import fred

def read_old_fredbin(filename):

    data = None

    # try reading in single line files
    try:
        dtype={'names': fred.fredbin_col_names,'formats': fred.fredbin_col_types}
        data = np.fromfile(filename, dtype)
        return data
    except:
        pass

    return None

def upgrade_fredbin(filename):

    # first try reading using standard methods
    data = fred.read_fredbin(filename)
    if data is not None:
        print("%s is up to date" % (filename))
        return None

    # read using other techniques
    data = read_old_fredbin(filename)

    if data is None:
        print("Cannot read %s " % (filename))
        return None

    # save the file
    with open(filename, 'wb') as outfile:
        print("Upgrading %s" % (filename))
        fred.write_fredbin(outfile, data)

    return filename


def main():
    parser = argparse.ArgumentParser(description='Upgrades fredbin files to the latest format.')
    parser.add_argument('input', nargs='+', help="FREDBIN files to be upgrade")
    parser.add_argument('-j','--jobs', type=int, help="Parallel jobs", default=4)
    parser.add_argument('--debug', default=False, action='store_true',
                        help="Run in debug mode")

    args = parser.parse_args()

    if args.debug:
        for filename in args.input:
            r = upgrade_fredbin(filename)
            results.extend(r)
    else:
        pool = mp.Pool(args.jobs)

        # farm out the jobs and wait for the result
        pool_result = pool.imap(upgrade_fredbin, args.input)
        pool.close()
        pool.join()

if __name__ == "__main__":
    main()
