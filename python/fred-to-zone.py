#!/bin/python
import argparse
import numpy as np
from numpy import cos
from math import pi
import inspect
import glob

# multithreading
from multiprocessing import Pool
from threading import Lock
import signal

# function parameter mapping
from functools import partial

# APASS-specific things
from quadtree import *
from quadtree_types import *
from apass import *
from apass_types import *

class WriteLock():
    """A class which provides mutually exclusive access to a file without using
    lockfiles"""

    def __init__(self, datapath, file_id):
        """"Create the WriteLock.

        datapath -- The full path to a directory where .dat files may be written
        file_id -- An integer representing a specific file"""
        self.filename = datapath + "/" + str(file_id).zfill(5) + ".dat"
        self.lock = Lock()
        self.outfile = None
        self.count = 0

    def open(self):
        """Opens the file or increment the reference counter"""
        with self.lock:

            if self.outfile is None:
                self.outfile = open(self.filename, 'a+b')

            self.count += 1

    def close(self):
        """Decrements the reference counter and closes the file when it is no
        longer needed"""
        with self.lock:
            self.count -= 1

            if self.count == 0:
                self.outfile.close()
                self.outfile = None

    def write(self, data):
        """Writes data to the file."""
        with self.lock:
            self.outfile.write(data)

class FileStore():
    """A class which manages multiple access to multiple files for multiple
    threads"""
    files = {}

    def __init__(self, datapath):
        """Create the Filestore

        datapath -- The full path to a directory to which files may be written."""
        self.datapath = datapath

    def open(self, file_id):
        """Open the file corresponding to file_id

        file_id -- an integer representing the file we wish to open"""
        if not file_id in self.files.keys():
            self.files[file_id] = WriteLock(self.datapath, file_id)

        self.files[file_id].open()

    def close(self, file_id):
        """Closes the file with the specified file_id"""
        self.files[file_id].close()

    def write(self, file_id, data):
        """Writes data to the file with the corresponding file_id"""
        self.files[file_id].write(data)


def fred_to_zone(tree_file, filestore, filename):
    """Processes an APASS FRED file into zones using data contained in the tree

    tree_file --- the full path to a JSON file containing tree generated using make-zones.py
    filestore --- a reference to a FileStore object
    filename --- the APASS FRED file to process"""
#    print("starting with %s %s" % (tree_file, filename))

    # restore the tree. make-zones.py writes out leaves of type IDLeaf
    tree = QuadTreeNode.from_file(tree_file, leafClass=IDLeaf)

    # read the fred file into a numpy array
    print("Processing FRED file " + filename)
    data = read_fred(filename)

    # process every file, inserting it into a .dat file. Keep track of any
    # files that were opened in the open_files set
    open_files = set()
    for datum in np.nditer(data):
        # pull out the RA and DEC
        ra  = datum['ra']
        dec = datum['dec']
        file_id = tree.insert(ra * cos(dec * pi / 180), dec, None)

        if file_id not in open_files:
            filestore.open(file_id)
            open_files.add(file_id)

        filestore.write(file_id, datum)

    # inform the filestore that we are done with the files.
    for file_id in open_files:
        filestore.close(file_id)

filestore = None

def main():

    parser = argparse.ArgumentParser(description='Merges .fred photometric files')
    #parser.add_argument('outdir', help="Directory into which .fredbin files will be generated")
    parser.add_argument('tree', help="Global zone map file produced by make-zones.py")
    parser.add_argument('input', nargs='+', help="Input files which will be split into zonefiles")
    parser.add_argument('-j','--jobs', type=int, help="Parallel jobs")
    parser.set_defaults(jobs=1)

    args = parser.parse_args()

    outdir = '/tmp/apass'
    wait_period = 120
    filestore = FileStore(outdir)

    # check for existing lockfiles
    lockfiles = glob.glob(outdir + "/*.lock")
    if len(lockfiles) > 0:
        print("Output directory contains lockfiles. Please remove to proceed")
        exit(0)

    # If you are going to profile the code, run this instead:
    #fred_to_zone(args.tree, filestore, args.input[0])
    #return

    # generate a pool of threads to process the input FRED files
    # dispatch to each thread using threadfunc
    pool = Pool(args.jobs)
    try:
        threadfunc = partial(fred_to_zone, args.tree, filestore)
        res = pool.map_async(threadfunc, args.input)
        res.get(wait_period) # wait 60 seconds for the thread to complete
    except KeyboardInterrupt:
        print("Caught keyboard interrupt, terminating workers.")
        pool.terminate()
    else:
        pool.close()

    pool.join()


if __name__ == "__main__":
    main()
