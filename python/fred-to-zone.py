#!/bin/python
import argparse
import numpy as np
from numpy import cos
from math import pi
import inspect
import glob
from enum import Enum

# parallel processing
import multiprocessing as mp

# APASS-specific things
from quadtree import *
from quadtree_types import *
from apass import *

class Operation(Enum):
    STOP = 0
    OPEN = 1
    CLOSE = 2
    WRITE = 3

class FileHandle:
    def __init__(self, zone_id):
        self.filename = apass_save_dir + "/" + name_zone_file(zone_id)
        self.count = 0
        self.outfile = None

    def open(self):
        if self.outfile is None:
            self.outfile = open(self.filename, 'a+b')
            self.count = 1
        else:
            self.count += 1

    def close(self):
        self.count -= 1
        if self.count <= 0:
            self.outfile.close()

    def write(self, data):
        self.outfile.write(data)

    def flush(self):
        self.outfile.flush()


def writer_func(queue):

    print("Starting writer process")
    handles = {}

    while(True):
        # pull an item off the queue
        operation, zone_id, data = queue.get()
        # print [operation, zone_id]

        # termination condition
        if operation == Operation.STOP:
            break

        if operation == Operation.OPEN:
            handle = FileHandle(zone_id)
            handle.open()
            handles[zone_id] = handle
        elif operation == Operation.CLOSE:
            handle = handles[zone_id]
            handle.flush()
            handle.close()
            if handle.outfile.closed:
                del handles[zone_id]
        elif operation == Operation.WRITE:
            handle = handles[zone_id]
            handle.write(data)
        else:
            raise RuntimeError("The requested operation is not implemented!")


def fred_to_zone_func(queue, filename):
    """Processes an APASS FRED file into zones using data contained in the tree

    tree_file --- the full path to a JSON file containing tree generated using make-zones.py
    filestore --- a reference to a FileStore object
    filename --- the APASS FRED file to process"""

    global tree_file

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
        [ra, dec] = get_coords(datum)
        zone_id = tree.insert(ra, dec, None)

        if zone_id not in open_files:
            enqueue_open(queue, zone_id)
            open_files.add(zone_id)

        enqueue_write(queue, zone_id, datum)

    # inform the filestore that we are done with the files.
    for zone_id in open_files:
        enqueue_close(queue, zone_id)

def enqueue_stop(queue):
    queue.put([Operation.STOP, -1, None])

def enqueue_open(queue, zone_id):
    queue.put([Operation.OPEN, zone_id, None])

def enqueue_close(queue, zone_id):
    queue.put([Operation.CLOSE, zone_id, None])

def enqueue_write(queue, zone_id, data):
    queue.put([Operation.WRITE, zone_id, data])


def main():

    parser = argparse.ArgumentParser(description='Merges .fred photometric files')
    #parser.add_argument('outdir', help="Directory into which .fredbin files will be generated")
    parser.add_argument('input', nargs='+', help="Input files which will be split into zonefiles")
    parser.add_argument('-j','--jobs', type=int, help="Parallel jobs")
    parser.set_defaults(jobs=1)

    args = parser.parse_args()
    outdir = apass_save_dir

    # load the global tree
    global tree_file
    tree_file = apass_save_dir + "/global.json"

    # set up the pool
    pool = mp.Pool(args.jobs)
    manager = mp.Manager()
    writer_queue = manager.Queue()

    # start up the writer
    writer_proc = mp.Process(target=writer_func, args=(writer_queue,))
    writer_proc.start()

    for filename in args.input:
        pool.apply_async(fred_to_zone_func, (writer_queue, filename))

    # wait for the pool to complete
    pool.close()
    pool.join()
    enqueue_stop()
    writer_proc.join()


if __name__ == "__main__":
    main()
