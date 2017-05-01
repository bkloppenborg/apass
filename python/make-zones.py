#!/bin/python
import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from numpy import cos
from math import pi

from quadtree import *
from quadtree_types import *

import apass

def merge_polar_zones(node):
    """Replaces QuadTree nodes that reside completely within the polar zones
    with a IDLeaf node with fileid = [1=North,2=South]"""

    north = 88
    south = -1 * north

    for i in range(len(node.children) - 1, -1, -1):
        child = node.children[i]
        rect = child.rect

        if rect.y_min == -90 and rect.y_max < south:
            rect = Rect(0, 360, -90, rect.y_max)
            node.children[i] = IDLeaf(rect, child.depth, node_id=2, parent=node)
        elif rect.y_max == 90 and rect.y_min > north:
            rect = Rect(0, 360, rect.y_min, 90)
            node.children[i] = IDLeaf(rect, child.depth, node_id=1, parent=node)

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
    fileid = 3 # reserve 1, 2 for the poles

    zonefile = apass.apass_save_dir + '/global.json'

    # subdivde the sphere to this depth:
    #depth = 6 # dRA = 5.625 dDEC = 2.8125
    depth = 7 # dRA = 2.8125 dDEC = 1.40625

    # build the quadtree, then merge the zones
    bounds = Rect(0, 360, -90, 90)
    tree = QuadTreeNode(bounds, 0)
    tree.split_until(depth, leafClass=IDLeaf)
    tree.runFunc(merge_polar_zones)

    if args.plot:
        plot_zones(tree)

    leaves = tree.get_leaves()
    print("Created %i zones" % (len(leaves)))

    # write the quadtree
    QuadTreeNode.to_file(tree, zonefile)

if __name__ == "__main__":
    main()
