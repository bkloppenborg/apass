#!/bin/python
import argparse
import numpy as np
import os
import time

#matplotlib
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches

# multithreading
from multiprocessing import Pool
from threading import Lock
import signal
from functools import partial

# APASS-specific things
from quadtree import *
from quadtree_types import *
from apass_types import *
import apass
from fred import read_fredbin
from border_info import make_border_info, save_border_info
import zone

def zone_to_rects(save_dir, filename):
    """Processes and APASS zone file into overlapping rectangles"""
    zone_id = apass.zone_from_name(filename)

    global tree_file

    # read in the (binary) data file
    data = read_fredbin(filename)
    print "Processing '%s' which has %i data points " % (filename, data.size)

    # find the bounds of this zone using the first data point in the file
    global_tree = QuadTreeNode.from_file(tree_file, leafClass=IDLeaf)
    datum = data[0]
    ra, dec = apass.get_coords(datum)
    zone_node = global_tree.find_leaf(ra, dec)
    zone_bounds = zone_node.rect

    # build a tree for the zone
    zone_tree = QuadTreeNode(zone_bounds, 0, parent=None)
    zone_tree.split_until(apass.zone_depth, leafClass=RectLeaf)

    # insert the data into the tree, building up containers (rectangles) in the
    # process
    for datum in np.nditer(data):
        datum = datum.copy()
        ra, dec = apass.get_coords(datum)
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
    zone_border_info = merge_containers_on_borders(zone_tree)
    #plot_rects(zone_tree) # plot the nodes after the merge

    # write out the containers that were on the border
    filename = filename_no_ext + '-border-rects.json'
    save_border_info(filename, zone_border_info)

    zone.save_zone_data(zone_tree, save_dir)

    # save the zone -> container mapping
    QuadTreeNode.to_file(zone_tree, filename_no_ext + "-zone.json")

def merge_containers_on_borders(zone_tree):

    zone_border_rects = dict()

    leaves = zone_tree.get_leaves() # these are of type RectLeaf
    for leaf in leaves:

        #  iterate over every RectContainer instance in the RectLeaf
        for container in leaf.containers:

            # get a unique list of nodes containing the corners
            # of this container. The polar zones return all other
            # polar zones so that we can collapse the data into
            # a single RectContainer instance
            adj_leaves = []
            if container.rect.y_min < -90:
                adj_leaves = get_polar_leaves(zone_tree, -90)
            elif container.rect.y_max > 90:
                adj_leaves = get_polar_leaves(zone_tree, 90)
            else:
                adj_leaves = get_containing_nodes(leaf, container)

            # remove self references
            adj_leaves.remove(leaf)

            for adj_leaf in adj_leaves:
                if adj_leaf is None:
                    # if there is no adjacent leaf, then we are at a zone border.
                    # make an entry in the zone_border_rects
                    info = make_border_info(container)
                    zone_border_rects.update(info)
                else:
                    # get the overlapping containers
                    adj_containers = adj_leaf.get_overlapping_containers(container)
                    # move their data
                    for adj_container in adj_containers:
                        container.merge(adj_container)
                        adj_leaf.remove_container(adj_container)

    return zone_border_rects

def get_polar_leaves(zone_tree, y_val):
    """Finds and returns RectLeaf nodes within the tree
    located at the poles.
    NOTE: Specify y_val as +90 or -90"""

    polar_leaves = []

    leaves = zone_tree.get_leaves()
    for leaf in leaves:
        if leaf.rect.y_min <= y_val or leaf.rect.y_max >= y_val:
            polar_leaves.append(leaf)

    return polar_leaves

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


def get_containing_nodes(leaf, container, ignore_self=False):
    """Return a list of nodes which contain the specified rectangle. The query
    originates from the specified leaf node and can move freely throughout the
    tree. This function may return None which should be interpreted as
    'the rectangle extends outside of the tree containing leaf'
    """
    nodes = []

    corners = container.get_corners()
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
    start = time.time()

    save_dir = os.path.dirname(os.path.realpath(args.input[0]))

    global tree_file
    tree_file = save_dir + "/global.json"

    # Construct a partial to serve as the function to call in serial or
    # parallel mode below.
    ztr_func = partial(zone_to_rects, save_dir)

    # use this for single thread development and debugging
    if args.debug:
        for filename in args.input:
            ztr_func(filename)
    else:
        # generate a pool of threads to process the input
        pool = Pool(args.jobs)

        # farm out the work
        result = pool.imap(ztr_func, args.input)
        pool.close()
        pool.join()

    # write out a file containing information on the containers modified.
    mod_file = save_dir + "/zone-to-rects-modified-files.txt"
    with open(mod_file, 'w') as outfile:
        for filename in args.input:
            path,fredbin_filename = os.path.split(filename)
            zone_id = apass.zone_from_name(fredbin_filename)
            filename = apass.name_zone_container_file(zone_id)
            outfile.write(save_dir + filename + "\n")

    print("A list of modified files has been written to %s" % (mod_file))

    end = time.time()
    print("Time elapsed: %is" % (int(end - start)))

if __name__ == "__main__":
    main()
