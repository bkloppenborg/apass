#!/bin/python
import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import os

# multithreading
from multiprocessing import Pool
from threading import Lock
import signal
from functools import partial

# APASS-specific things
from quadtree import *
from quadtree_types import *
from apass import *
from apass_types import *


def zone_to_rects(filename):
    """Processes and APASS zone file into overlapping rectangles"""
    zone_id = zone_from_name(filename)

    global tree_file
    depth = 4

    # read in the (binary) data file
    data = read_fredbin(filename)
    print "Processing '%s' which has %i data points " % (filename, data.size)

    # find the bounds of this zone using the first data point in the file
    global_tree = QuadTreeNode.from_file(tree_file, leafClass=IDLeaf)
    datum = data[0]
    ra, dec = get_coords(datum)
    zone_node = global_tree.find_leaf(ra, dec)
    zone_bounds = zone_node.rect

    # build a tree for the zone
    zone_tree = QuadTreeNode(zone_bounds, 0, parent=None)
    zone_tree.split_until(depth, leafClass=RectLeaf)

    # insert the data into the tree, building up containers (rectangles) in the
    # process
    for datum in np.nditer(data):
        datum = datum.copy()
        ra, dec = get_coords(datum)
        try:
            zone_tree.insert(ra, dec, datum)
        except RuntimeError:
            print("ERROR: Potential data corruption in " + filename)
            print("ERROR: Check file, remove the zone directory, and re-run this program")
            return


    # prepare the save the data.
    # the zone file's name
    filename_no_ext = os.path.splitext(filename)[0]

    # now number the (leaf) nodes
    number_containers(zone_tree, zone_id=zone_id)
    #plot_rects(zone_tree) # plot the nodes before the merge
    leaves = zone_tree.get_leaves()
    zone_border_info = merge_containers_on_borders(leaves)
    #plot_rects(zone_tree) # plot the nodes after the merge

    # write out the containers that were on the border
    filename = filename_no_ext + '-border-rects.json'
    save_border_info(filename, zone_border_info)

    apass.save_zone_data(zone_tree, apass.apass_save_dir)

    # save the zone -> container mapping
    QuadTreeNode.to_file(zone_tree, filename_no_ext + "-zone.json")

def number_containers(tree, zone_id=0):
    """Assigns a unique ID to each container."""

    node_id = 0
    leaves = tree.get_leaves()

    # number the leaves with unique node identifiers
    for leaf in leaves:
        if not isinstance(leaf, RectLeaf):
            raise RuntimeError("Encountered an unexpected non-RectLeaf node!")

        leaf.zone_id = zone_id
        leaf.node_id = node_id

        # number the containers with unique container IDs
        container_id = 0
        for container in leaf.containers:
            container.zone_id = zone_id
            container.node_id = node_id
            container.container_id = container_id

            # number the data with the corresponding (node_id, container_id)
            for datum in container.data:
                datum['node_id'] = node_id
                datum['container_id'] = container_id

            # increment the container ID
            container_id += 1

        # increment the node ID
        node_id += 1

def merge_containers_on_borders(nodes):
    """This function merges containers that reside on the border of cells.
    Any container that resides on the border of the zone will be returned
    as a border_info dict made by apass.make_border_info(...)"""

    zone_border_rects = dict()
    for node in nodes:

        if not isinstance(node, RectLeaf):
            continue

        # iterate over every RectContainer in the RectLeaf node
        for container in node.containers:
            rect = container.rect

            # Get a unique list of nodes containing the corners of the rectangle
            # excluding this node. These will be of the type RectLeaf.
            other_nodes = get_containing_nodes(node, rect)
            other_nodes.remove(node)

            # iterate through the list of other nodes and attempt to merge
            # containers
            for other_node in other_nodes:
                if other_node is None:
                    info = apass.make_border_info(container)
                    zone_border_rects.update(info)
                else:
                    other_containers = other_node.get_overlapping_containers(container, remove=True)
                    for other in other_containers:
                        container.merge(other)

    return zone_border_rects

def get_containing_nodes(leaf, rect, ignore_self=False):
    """Return a list of nodes which contain the specified rectangle. The query
    originates from the specified leaf node and can move freely throughout the
    tree. This function may return None which should be interpreted as
    'the rectangle extends outside of the tree containing leaf'
    """
    nodes = []
    corners = rect.get_corners()
    for x,y in corners:

        if leaf.contains(x,y):
            if leaf not in nodes:
                nodes.append(leaf)
        else:
            # This function call can return None, which is a valid response.
            # This indicates that the containing node is not in this tree.
            other = leaf.find_node_containing(x,y)
            if other not in nodes:
                nodes.append(other)

    return nodes

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
    parser.add_argument('-j','--jobs', type=int, help="Parallel jobs", default=4)
    parser.add_argument('--debug', default=False, action='store_true',
                        help="Run in debug mode")
    parser.set_defaults(jobs=1)

    args = parser.parse_args()

    global tree_file
    tree_file = apass_save_dir + "/global.json"

    # use this for single thread development and debugging
    if args.debug:
        for filename in args.input:
            zone_to_rects(filename)
    else:
        # generate a pool of threads to process the input
        pool = Pool(args.jobs)

        # farm out the work
        result = pool.imap(zone_to_rects, args.input)
        pool.close()
        pool.join()

if __name__ == "__main__":
    main()
