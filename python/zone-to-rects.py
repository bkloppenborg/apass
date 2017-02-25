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

    # find the bounds of this zone using the first data point in the file
    global_tree = QuadTreeNode.from_file(tree_file, leafClass=IDLeaf)
    datum = data[0]
    ra, dec = get_coords(datum)
    zone_node = global_tree.find_leaf(ra, dec)
    zone_bounds = zone_node.rect
    #print zone_bounds

    # build a tree for the zone
    print("Building tree")
    zone_tree = QuadTreeNode(zone_bounds, 0, parent=None)
    zone_tree.split_until(depth, leafClass=RectLeaf)

    # insert the data into the tree, building up containers (rectangles) in the
    # process
    print("Inserting data into the tree")
    print("Data size: %i" % (data.size))
    for datum in np.nditer(data):
        ra, dec = get_coords(datum)
        zone_tree.insert(ra, dec, datum)


    # prepare the save the data. Begin by setting the output directory to match
    # the zone file's name
    directory = os.path.splitext(filename)[0]
    if not os.path.exists(directory):
        os.makedirs(directory)

    # now number the (leaf) nodes
    number_containers(zone_tree, zone_id)
    #plot_rects(zone_tree) # plot the nodes before the merge
    zone_border_containers = zone_tree.runFuncRet(merge_containers_on_borders)
    #plot_rects(zone_tree) # plot the nodes after the merge

    # write out the containers that were on the border
    filename = directory + '/border-rects.txt'
    outfile = open(filename, 'w')
    for container in zone_border_containers:
        name = apass.name_rect(container.zone_id, container.node_id, container.container_id)
        outfile.write(name + "\n")
    outfile.close()

    # save the rectangle data to file
    func = partial(save_data, directory=directory)
    zone_tree.runFunc(func)

    # save the zone -> container mapping
    QuadTreeNode.to_file(zone_tree, directory + "/zone.json")

def save_data(node, directory="/tmp"):
    """Saves the data from the RectLeaf's RectContainers to disk."""

    # this runs only on RectLeaf objects
    if not isinstance(node, RectLeaf):
        return

    # write out the data in each container.
    for container in node.containers:
        container.save(directory)

    # remove the RectContainers from this node for serialization
    node.containers = []

def number_containers(tree, zone_id=0):

    node_id = 0
    leaves = tree.get_leaves()

    for leaf in leaves:
        if not isinstance(leaf, RectLeaf):
            raise RuntimeError("Encountered an unexpected non-RectLeaf node!")

        leaf.node_id = node_id;
        node_id += 1

        container_id = 0
        for container in leaf.containers:
            container.zone_id = zone_id
            container.node_id = node_id
            container.container_id = container_id
            container_id += 1

def merge_containers_on_borders(node):
    """This function merges containers that reside on the border of cells.
    Any container that resides on the border of the zone will be returned
    from this function in the format [zone_id, node_id, container_id]"""

    zone_border_rects = []

    if not isinstance(node, RectLeaf):
        return zone_border_rects

    # iterate over every RectContainer in the RectLeaf node
    for container in node.containers:
        rect = container.rect

        # Get a unique list of nodes containing the corners of the rectangle.
        # These will be of the type RectLeaf
        nodes = get_containing_nodes(node, rect)
        nodes.remove(node)

        # iterate through the list of other nodes and attempt to merge
        # containers
        for other_node in nodes:
            # If a node is None, that means the rectangle is on the border of a
            # zone. Save this RectContainer for later analysis
            if other_node is None:
                if node not in zone_border_rects:
                    zone_border_rects.append(container)
            else:
                other_node.merge_containers(container)

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
    parser.add_argument('-j','--jobs', type=int, help="Parallel jobs")
    parser.set_defaults(jobs=1)

    args = parser.parse_args()

    global tree_file
    tree_file = apass_save_dir + "/global.json"

    # use this for single thread development and debugging
    #for filename in args.input:
    #    zone_to_rects(filename)
    #return

    # generate a pool of threads to process the input FRED files
    # dispatch to each thread using threadfunc
    pool = Pool(args.jobs)
    try:
        # use the following during production runs, keyboard interrupts won't be handled
        pool.map(zone_to_rects, args.input)
        # use this when developing/debugging, interrupts handled
        #res = pool.map_async(zone_to_rects, args.input)
        #res.get(wait_period)
    except KeyboardInterrupt:
        print("Caught keyboard interrupt, terminating workers.")
        pool.terminate()
    else:
        pool.close()

    pool.join()

if __name__ == "__main__":
    main()
