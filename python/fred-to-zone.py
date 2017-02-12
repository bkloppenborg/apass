#!/bin/python
import argparse
import numpy as np
from numpy import cos
from math import pi

from functools import partial
import inspect

from quadtree import *
from quadtree_types import *
from apass import *
from apass_types import *

def set_filestores(node, filestores=None):
    if isinstance(node, FileStoreLeaf):
        node.filestores = filestores

def main():

    parser = argparse.ArgumentParser(description='Merges .fred photometric files')
    #parser.add_argument('outdir', help="Directory into which .fredbin files will be generated")
    parser.add_argument('tree', nargs=1, help="Global zone map file produced by make-zones.py")
    parser.add_argument('input', nargs='+', help="Input files which will be split into zonefiles")

    args = parser.parse_args()

    filestores = FileStores("/tmp/apass/")

    # Load the quadtree, this file must exist
    tree = None
    #try:
    tree = QuadTreeNode.from_file('/tmp/apass/global.json', leafClass=FileStoreLeaf)
    partial_func = partial(set_filestores, filestores=filestores)
    tree.runFunc(partial_func)
    #except:
    #    print("Could not load global quadtree! You must create this file first with make-zones.py")
    #    return 1

    for filename in args.input:
        print("Reading FRED file " + filename)
        data = read_fred(filename)
        print("Finished Reading FRED")

        for datum in np.nditer(data):
           # pull out the RA and DEC
           ra  = datum['ra']
           dec = datum['dec']
           tree.insert(ra * cos(dec * pi / 180), dec, datum)

        filestores.close()

if __name__ == "__main__":
    main()
