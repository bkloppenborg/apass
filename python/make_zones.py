#!/bin/python
import argparse
import numpy as np
import time

# matplotlib
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches

# various mathematical functions
from numpy import cos
from math import pi

from quadtree import *
from quadtree_types import *

import apass

def merge_polar_zones(root_node):
    """Replaces QuadTree nodes that reside completely within the polar zones
    with a IDLeaf node with fileid = [1=North,2=South]"""

    north = apass.polar_zone_cutoff
    south = -1 * north

    leaves = root_node.get_leaves()

    for leaf in leaves:
        rect = leaf.rect

        if rect.y_min == -90 or rect.y_max < south:
            leaf.node_id = apass.south_zone_id
            leaf.rect.x_min = 0
            leaf.rect.x_max = 360
        elif rect.y_max == 90 or rect.y_min > north:
            leaf.node_id = apass.north_zone_id
            leaf.rect.x_min = 0
            leaf.rect.x_max = 360

def number_zones(root_node):
    """Assigns a unique ID to each zone following a raster-scan pattern"""

    counter = max(apass.north_zone_id, apass.south_zone_id) + 1

    x_max = 2**apass.global_depth
    y_max = 2**apass.global_depth

    dx = 360.0 / x_max
    dy = 180.0 / y_max

    for j in range(0, y_max):
        for i in range(0, x_max):

            x = dx * i
            y = dy * j - 90.0

            node = root_node.find_leaf(x,y)

            if node.node_id == None:
                node.node_id = counter
                counter += 1

    leaves = root_node.get_leaves()
    for leaf in leaves:
        leaf.node_id = counter
        counter += 1


def export_rect(node):
    """Export leaf rectangles to the global rects variable"""
    if node.is_leaf():
        global rects
        rects.append(node.rect)

def plot_zones(save_dir, tree):
    """Plots leaf zones found in the global rects variable"""

    fig = plt.figure(figsize=(120,100))
    axes = plt.gca()
    axes.set_xlim([0,360])
    axes.set_ylim([-90,90])

    leaves = tree.get_leaves()
    for leaf in leaves:
        # extract the bounds of this leaf
        x_min = leaf.rect.x_min
        x_max = leaf.rect.x_max
        y_min = leaf.rect.y_min
        y_max = leaf.rect.y_max

        # calculate the width, height, and center
        width  = x_max - x_min
        height = y_max -  y_min
        center_x = (x_max + x_min) / 2
        center_y = (y_max + y_min) / 2

        # plot a rectangle for this zone and number it
        axes.add_patch(patches.Rectangle((x_min,y_min), width, height, fill=False))
        axes.text(center_x, center_y, leaf.node_id)

    plt.savefig(save_dir + "/global.png")

rects = [] # stores exported rectangles

def main():

    parser = argparse.ArgumentParser(description='Merges .fred photometric files')
    parser.add_argument('save_dir', help="Directory to save the output files.")
    parser.add_argument('--plot', dest='plot', help="Plot the generated zones", action='store_true')
    parser.set_defaults(plot=False)

    # parse the command line arguments and start timing the script
    args = parser.parse_args()
    start = time.time()

    global fileid # used in quadtree_types.py
    fileid = apass.north_zone_id + 1 # skip over reserved IDs

    zonefile = args.save_dir + '/global.json'

    # build the quadtree, then merge the zones
    bounds = Rect(0, 360, -90, 90)
    tree = QuadTreeNode(bounds, 0)
    tree.split_until(apass.global_depth, leafClass=IDLeaf)
    number_zones(tree)
    merge_polar_zones(tree)

    if args.plot:
        plot_zones(args.save_dir, tree)

    leaves = tree.get_leaves()
    print("Created %i zones" % (len(leaves)))

    # write the quadtree
    QuadTreeNode.to_file(tree, zonefile)

    end = time.time()
    print("Time elapsed: %is" % (int(end - start)))

if __name__ == "__main__":
    main()
