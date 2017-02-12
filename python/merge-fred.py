#!/bin/python
import sys
import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from numpy import cos
from math import pi

# Custom modules
sys.path.append(os.path.abspath('./modules/FileLock/filelock'))
from filelock import FileLock

from quadtree import *
from apass import *

class FileStore():

    def __init__(self, filename):
        self.filename = filename

    def close(self):
        self.file.close()
        self.lock.release()

    def open(self):
        # obtain a lock or keep retrying for 60 seconds
        self.lock = FileLock(self.filename, timeout=60, delay=0.05)
        self.lock.acquire()
        self.file = open(self.filename, 'a+')

    def append(self, data):
        self.file.write(data)

class FileStores():
    files = {}

    def __init__(self, datapath):
        self.datapath = datapath

    def close(self):
        for id, filestore in self.files.iteritems():
            filestore.close()

    def get_filehandle(self, file_id):

        if file_id not in self.files:
            filename = self.datapath + "/" + str(file_id).zfill(5) + ".dat"
            filestore = FileStore(filename)
            filestore.open()
            self.files[file_id] = filestore

        return self.files[file_id]

    def insert(self, file_id, data):
        filehandle = self.get_filehandle(file_id)
        filehandle.append(data)


class FileStoreLeaf(QuadTreeNode):
    """Class for a QuadTree leaf that points to a file"""

    def __init__(self, rect, depth, file_id=None, parent=None):
        QuadTreeNode.__init__(self, rect, depth, parent)

        self.file_id = file_id;
        if self.file_id is None:
            global fileid
            self.file_id = fileid
            fileid += 1

    @staticmethod
    def from_dict(rect, depth, dict_):
        file_id = dict_['file_id']
        return FileStoreLeaf(rect, depth, file_id=file_id)

    def insert(self, x, y, data):
        """Inserts the data into the corresponding file"""
        global filestores
        filestores.insert(self.file_id, data)


def merge_polar_zones(node):
    """Replaces QuadTree nodes that reside completely within the polar zones
    with a FileStoreLeaf node with fileid = [0=North,1=South]"""

    north = 80
    south = -1 * north

    for i in range(len(node.children) - 1, -1, -1):
        child = node.children[i]
        rect = child.rect

        if rect.y_min == -90 and rect.y_max < south:
            rect = Rect(0, 360, -90, rect.y_max)
            node.children[i] = FileStoreLeaf(rect, child.depth, file_id=1, parent=node)
        elif rect.y_max == 90 and rect.y_min > north:
            rect = Rect(0, 360, rect.y_min, 90)
            node.children[i] = FileStoreLeaf(rect, child.depth, file_id=0, parent=node)

def export_rect(node):
    if node.is_leaf():
        global rects
        rects.append(node.rect)

def plot_zones(tree):
    tree.runFunc(export_rect)

    fig = plt.figure()
    axes = plt.gca()
    axes.set_xlim([0, 360])
    axes.set_ylim([-90, 90])
    for rect in rects:
        x = rect.x_min
        y = rect.y_min
        width = rect.x_max - rect.x_min
        height = rect.y_max -  rect.y_min
        axes.add_patch(patches.Rectangle((x,y), width, height, fill=False))

    plt.show()


filestores = None
fileid = 2 # reserve 0 and 1 for the pole
rects = [] # stores exported rectangles

def main():

    parser = argparse.ArgumentParser(description='Merges .fred photometric files')
    #parser.add_argument('outdir', help="Directory into which .fredbin files will be generated")
    parser.add_argument('input', nargs='+', help="Input files which will be merged")

    args = parser.parse_args()

    global filestores
    filestores = FileStores("/tmp/apass/")

    # comes from command line
    depth = 6 # dRA = 5.625 dDEC = 2.8125
    #depth = 7 # dRA = 2.8125 dDEC = 1.40625

    bounds = Rect(0, 360, -90, 90)
    tree = QuadTreeNode(bounds, 0)
    tree.split_until(depth, leafClass=FileStoreLeaf)
    tree.runFunc(merge_polar_zones)

    QuadTreeNode.to_file(tree, '/tmp/apass/global.json')

    temp = QuadTreeNode.from_file('/tmp/apass/global.json', leafClass=FileStoreLeaf)

    print tree.size()
    print temp.size()

    quit()


    dtype={'names': apass_col_names,'formats': apass_col_types}
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
