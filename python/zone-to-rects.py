#!/bin/python
import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

# multithreading
#from multiprocessing import Pool
#from threading import Lock
#import signal

# APASS-specific things
from quadtree import *
from quadtree_types import *
from apass import *

class RectContainer():

    def __init__(self, x, y, data):

        dr = 1. / (60 * 60) # 1 arcsecond in degree
        dx = dr * cos(y * pi /  180)
        dy = dr
        self.rect = Rect(x - dx, x + dx, y - dy, y + dy)
        self.data = [data]

    def overlaps(self, other):
        """Determines if this RectContainer overlaps with another RectContainer

        other -- another RectContainer instance"""
        return self.rect.overlaps(other.rect)

    def merge(self, other):
        """Merges two RectContainer Instances, growing their bounding rectangles

        other -- another RectContainer instance"""
        self.rect.expand(other.rect)
        self.data.extend(other.data)

class RectLeaf(QuadTreeNode):
    """A class for a QuadTree leaf that contains a list of rectangles with
    data."""

    def __init__(self, rect, depth, parent=None):
        QuadTreeNode.__init__(self, rect, depth, parent)

        self.containers = []

    def insert(self, x, y, data):
        """Stores the data inside of a container encapsulated by this node."""

        # build a rectangle container for this data
        rc = RectContainer(x, y, data)

        # find, merge, and remove containers that overlap with rc
        containers = []
        for container in self.containers:
            if rc.overlaps(container):
                rc.merge(container)
            else:
                containers.append(container)

        containers.append(rc)
        self.containers = containers

def zone_to_rect(filename):
    """Processes and APASS zone file into overlapping rectangles"""
    global tree_file

    depth = 4

    # read in the (binary) data file
    data = read_fredbin(filename)

    # find the bounds of this zone using the first data point in the file
    global_tree = QuadTreeNode.from_file(tree_file, leafClass=IDLeaf)
    datum = data[0]
    ra, dec = get_coords(datum)
    zone_node = global_tree.find_leaf(ra, dec)
    zone_bounds = zone_node.rect
    print zone_bounds

    # build the zone tree
    print("Building tree")
    zone_tree = QuadTreeNode(zone_bounds, 0, parent=None)
    zone_tree.split_until(depth, leafClass=RectLeaf)

    # insert the data into the tree, building up a rectangle in the process
    print("Inserting data into the tree")
    print("Data size: %i" % (data.size))
    for datum in np.nditer(data):
        ra, dec = get_coords(datum)
        zone_tree.insert(ra, dec, datum)

    plot_rects(zone_tree)

def export_container_data(node):
    """Export leaf rectangles to the global rects variable."""
    if node.is_leaf():
        global cells
        global containers
        global data

        cells.append(node.rect)

        for container in node.containers:
            containers.append(container.rect)
            data.extend(container.data)

def plot_rects(tree):
    """Plots leaf zones found in the global rects variable"""
    global cells
    global containers
    global data
    cells = []
    containers = []
    data = []

    tree.runFunc(export_container_data)

    print("leaf cells: %i containers: %i data: %i" % (len(cells), len(containers), len(data)))

    bounds = tree.rect
    xlim = [bounds.x_min, bounds.x_max]
    ylim = [bounds.y_min, bounds.y_max]

    fig = plt.figure()
    axes = plt.gca()
    axes.set_xlim(xlim)
    axes.set_ylim(ylim)

    for rect in cells:
        x = rect.x_min
        y = rect.y_min
        width = rect.x_max - rect.x_min
        height = rect.y_max -  rect.y_min
        axes.add_patch(patches.Rectangle((x,y), width, height, fill=False, edgecolor="red"))

    for rect in containers:
        x = rect.x_min
        y = rect.y_min
        width = rect.x_max - rect.x_min
        height = rect.y_max -  rect.y_min
        axes.add_patch(patches.Rectangle((x,y), width, height, fill=False))

    for datum in data:
        ra, dec = get_coords(datum)
        plt.scatter(ra, dec)

    plt.show()

def main():

    parser = argparse.ArgumentParser(description='Merges .fred photometric files')
    parser.add_argument('input', nargs='+', help="Input files which will be split into zonefiles")
    parser.add_argument('-j','--jobs', type=int, help="Parallel jobs")
    parser.set_defaults(jobs=1)

    args = parser.parse_args()

    global tree_file
    tree_file = apass_save_dir + "/global.json"

    zone_to_rect(args.input[0])

if __name__ == "__main__":
    main()
