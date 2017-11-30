#!/bin/python
import argparse
import numpy as np

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
        elif rect.y_max == 90 or rect.y_min > north:
            leaf.node_id = apass.north_zone_id

def number_zones(root_node):

    counter = max(apass.north_zone_id, apass.south_zone_id) + 1

    leaves = root_node.get_leaves()
    for leaf in leaves:
        leaf.node_id = counter
        counter += 1


def export_rect(node):
    """Export leaf rectangles to the global rects variable"""
    if node.is_leaf():
        global rects
        rects.append(node.rect)

def plot_zones(tree):
    """Plots leaf zones found in the global rects variable"""
    tree.runFunc(export_rect)
    bounds = tree.rect
    xlim = [bounds.x_min, bounds.x_max]
    ylim = [bounds.y_min, bounds.y_max]

    fig = plt.figure()
    axes = plt.gca()
    axes.set_xlim(xlim)
    axes.set_ylim(ylim)

    for rect in rects:
        x = rect.x_min
        y = rect.y_min
        width = rect.x_max - rect.x_min
        height = rect.y_max -  rect.y_min
        axes.add_patch(patches.Rectangle((x,y), width, height, fill=False))

    plt.show()

rects = [] # stores exported rectangles

def main():

    parser = argparse.ArgumentParser(description='Merges .fred photometric files')
    #parser.add_argument('outdir', help="Directory into which .fredbin files will be generated")
    parser.add_argument('--plot', dest='plot', help="Plot the generated zones", action='store_true')
    parser.set_defaults(plot=False)

    args = parser.parse_args()

    global fileid # used in quadtree_types.py
    fileid = apass.north_zone_id + 1 # skip over reserved IDs

    zonefile = apass.apass_save_dir + '/global.json'

    # build the quadtree, then merge the zones
    bounds = Rect(0, 360, -90, 90)
    tree = QuadTreeNode(bounds, 0)
    tree.split_until(apass.global_depth, leafClass=IDLeaf)
    number_zones(tree)
    merge_polar_zones(tree)

    if args.plot:
        plot_zones(tree)

    leaves = tree.get_leaves()
    print("Created %i zones" % (len(leaves)))

    # write the quadtree
    QuadTreeNode.to_file(tree, zonefile)

if __name__ == "__main__":
    main()
