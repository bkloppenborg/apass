#!/usr/bin/python

# system includes
import argparse
import sys, os
from multiprocessing import Pool
from functools import partial
import time

# for plots
import numpy as np

# matplotlib
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

# local includes
import apass
from quadtree import QuadTreeNode
from quadtree_types import *
from zone import load_zone, save_zone

def get_active_indices(i, j, stride):
    """Returns list of indices selected from (i,j) in steps of stride in each direction."""

    global width
    global height

    i_vals = range(i, width, stride)
    j_vals = range(j, height, stride)

    output = []
    for i in i_vals:
        for j in j_vals:
            if j < 2 or j > height - 2:
                continue

            output.append((i,j))

    return output


def get_adjacent_zone_ids(i, j):
    """Returns zones that should be regarded as adjacent to the specified zone.

    This function behaves as follows:
     * A polar zone: Returns the entire adjacent row
     * First row before/after a polar zone: Six adjacent zones that are not the polar zone
     * All others: The eight zones surrounding the zone of interest
    """

    global global_tree
    global height
    global width

    #print("Finding adjacent nodes for %i %i" % (i, j))

    adj_zone_indices = []

    # North Polar Zone
    if j == 0:
        for k in range(0, width):
            adj_zone_indices.append((k, j + 1))

    # Southern Polar Zone
    elif j == height - 1:
        for k in range(0, width):
            adj_zone_indices.append((k, j - 1))

    # first row after a polar zone
    elif j == 1:
        for k in range(-1, 2):
            for l in range(0, 2):
                if k == 0 and l == 0:
                    continue
                adj_zone_indices.append((i + k, j + l))

    # first row before a polar zone
    elif j == height - 2:
        for k in range(-1, 2):
            for l in range(-1, 1):
                if k == 0 and l == 0:
                    continue
                adj_zone_indices.append((i + k, j + l))

    # all other cases
    # return a the 3x3 block surrounding (i,j)
    else:
        for k in range(-1, 2):
            for l in range(-1, 2):
                if k == 0 and l == 0:
                    continue

                adj_zone_indices.append((i + k, j + l))

    #print "Adjacent zones: "
    #print adj_zone_indices

    adj_zone_ids = []
    for k,l in adj_zone_indices:
        adj_zone_id = get_zone_id_from_indices(k, l)
        adj_zone_ids.append(adj_zone_id)
    return adj_zone_ids

def get_zone_id_from_indices(i,j):
    """Returns a reference to the zone in the global tree given
    the (i,j) indices of the zone in the tree.

    Requires that the width, height, and global_tree are loaded into
    the global namespace.
    """
    global width
    global height

    i = i % width
    j = j % height

    # compute the center of the zone
    dRA = float(360) / width
    dDEC = float(180) / height
    x = float(i) / width * 360 + dRA / 2
    y = float(j) / height * 180 - 90 + dDEC / 2

    #print i, j, x, y

    zone_id = get_zone_id_from_coordinates(x,y)
    return zone_id

def get_zone_id_from_coordinates(x,y):
    """Returns a reference to a zone in the global tree given the specified
    (x,y) coordinates.

    Requires the global_tree is loaded into the global namespace."""
    global global_tree
    leaf = global_tree.find_leaf(x,y)
    return leaf.node_id

def lookup_zone_ids(save_dir, zone_position):

    # Get this zone and its ID
    i,j = zone_position
    zone_id      = get_zone_id_from_indices(i,j)
    adj_zone_ids = get_adjacent_zone_ids(i, j)

    return zone_id, adj_zone_ids

def fix_overlaps_from_positions(save_dir, zone_position):
    """Wrapper to fix_overlaps when only the (i,j) position of the zone is known.

    save_dir - the directory containing zone files
    zone_position - a (i,j) integer tuple indicating the location of the primary
                    zone within the over-all problem space
    """

    # load the global tree
    global global_tree
    global_tree_filename = save_dir + "/global.json"
    global_tree = QuadTreeNode.from_file(global_tree_filename, leafClass=IDLeaf)

    zone_id, adjacent_zone_ids = lookup_zone_ids(save_dir, zone_position)

    fix_overlaps(save_dir, zone_id, adjacent_zone_ids)

def fix_overlaps(save_dir, zone_id, adjacent_zone_ids):
    """Fixes overlapping containers between adjacent zones.
    NOTE: Requires global_tree_filename to be loaded into the global namespace!

    save_dir - the directory containing zone files
    zone_id - ID (int) of primary zone of interest for overlap detection
    adjacent_zone_ids - IDs (list[int]) of adjacent zones
    """

    print("Processing Zone: %i" % (zone_id))

    global global_tree
    global_tree_filename = save_dir + "/global.json"
    global_tree = QuadTreeNode.from_file(global_tree_filename, leafClass=IDLeaf)

    # Load the zone
    zone_dict = load_zone(save_dir, zone_id)
    tree = zone_dict['tree']
    border_infos = zone_dict['border_info']

    if not zone_dict['loaded']:
        #print("Failed to load zone %i " % (zone_id))
        return
    if len(border_infos) == 0:
        #print("No overlapping data found in zone %i" % (zone_id))
        return

    # (attempt to) load adjacent zone data.
    # If this fails adj_zones[id]['loaded'] will be false
    adj_zone_dicts = dict()
    for adj_zone_id in adjacent_zone_ids:
        adj_zone_dicts[adj_zone_id] = load_zone(save_dir, adj_zone_id)

    # iterate over this zone's border rectangles and merge overlapping entries
    for key, info in border_infos.iteritems():
        #print(" checking for overlaps with %s " % (key))

        # Find this zone's node and container which store this border rectangle
        x,y = info['center']
        node = tree.find_leaf(x,y)
        node_id = node.node_id
        container = node.get_container(x,y)

        # produce a nice name for this container
        dest_name = apass.name_container(container.zone_id,
                                         container.node_id,
                                         container.container_id)

        # Get the coordinates for the corners of this container
        corners = container.get_corners()

        # find adjacent nodes which overlap with this container
        for c_x, c_y in corners:
            c_x, c_y = apass.wrap_bounds(c_x, c_y)

            adj_zone = global_tree.find_leaf(c_x, c_y)
            adj_zone_id = adj_zone.node_id

            # skip corners that are within this zone
            if adj_zone_id == zone_id:
                continue

            # skip zones that didn't load
            if adj_zone_id not in adj_zone_dicts:
                continue
            if adj_zone_dicts[adj_zone_id]['loaded'] == False:
                continue

            # pull out some variables for easier referencing
            adj_tree = adj_zone_dicts[adj_zone_id]['tree']
            adj_border_infos = adj_zone_dicts[adj_zone_id]['border_info']

            # find the adjacent node and any adjacent containers
            adj_node = adj_tree.find_leaf(c_x, c_y)
            adj_containers = adj_node.get_overlapping_containers(container)

            # move the adjacent container's data into this container
            for adj_container in adj_containers:
                # generate nice output
                src_name = apass.name_container(adj_container.zone_id,
                                            adj_container.node_id,
                                            adj_container.container_id)
                print(" merging %s into %s" % (src_name, dest_name))
                container.merge(adj_container)
                adj_node.remove_container(adj_container)

                # remove the adjacent container from the border info file
                if src_name in adj_border_infos.keys():
                    #print(" updating border rect file for zone %i" % (adj_zone_id))
                    del adj_border_infos[src_name]

    # All done with the merge. Write out data for all files.
    save_zone(save_dir, zone_dict)
    for key, value in adj_zone_dicts.iteritems():
        save_zone(save_dir, value)


def main():
    parser = argparse.ArgumentParser(
        description='Resolves instances where star data spans a zone boundary')
    parser.add_argument('save_dir', help="Directory to save the output files.")
    parser.add_argument('-j','--jobs', type=int, help="Parallel jobs", default=4)
    parser.add_argument('--debug', default=False, action='store_true',
                        help="Run in debug mode")
    parser.add_argument('--zone', type=int, nargs='+',
                        help="Zone IDs for primary and ajacent zones to inspect.")
    parser.set_defaults(jobs=1)

    # parse the command line arguments and start timing the script
    args = parser.parse_args()
    start = time.time()


    # determine the size of the image
    global width
    global height
    width = 2**(apass.global_depth)
    height = width

    # move data that sticks out past the left to the right
    # TODO: Run this in parallel
    move_zero_edge_data(args.save_dir)

    # If requested, process a single zone
    if args.zone:
        if len(args.zone) == 1:
            print("You must specify a zone and at least one adjacent zone")
            return

        zone_id = args.zone[0]
        adjacent_zone_ids = args.zone[1:]
        fix_overlaps(args.save_dir, zone_id, adjacent_zone_ids)
        return

    print("Image size: %i %i" % (width, height))

    stride = 3 # need to lock a 3x3 block surrounding the zone.

    # Processing below proceeds as follows:
    #  1. Both polar zones in parallel
    #  2. First/Last row of the zone grid in parallel.
    #  3. All remaining zones in parallel.
    # Each of these methods steps through and locks the zones differently
    # as is explained below.

    # 1. Both polar zones in parallel
    # This will lock the polar zone and ALL zones in the adjacent row. E.g.
    # |           p           |
    # |x|x|x|x|x|x|x|x|x|x|x|x|
    # where 'p' is the polar zone and 'o' denotes locked zones
    zone_indices = [(0,0), (0, height - 1)]
    print("\nProcessing polar zones \n")
    process_indices(zone_indices, args)

    # 2. First/Last row of the zone grid in parallel.
    # This steps through the first and last non-polar row in strides of 3
    # (not counting the zeroth entry). E.g.:
    # |           p           |
    # |o|x|x|o|x|x|o|x|...|x|o|
    # |x|x|x|x|x|x|x|x|...|x|x|
    # where 'o' are the zones executed in parallel and 'x' denote locked zones
    # and 'p' is the polar zone. The 'o' entries are advanced one entry to
    # the right each time the for loop is advanced
    print("\nProcessing first/last row\n")
    for i in range(0, 3):
        zone_indices = []

        cols = range(i, width, stride)
        for col in cols:
            zone_indices.append((col, 1))

        for col in cols:
            zone_indices.append((col, height-2))

        process_indices(zone_indices, args)

    # 3. All remaining zones in parallel
    # Here we process all non-polar non-first/last row zones in parallel in
    # strides of 3 in row and 3 in column. The zeroth invocation will look like
    # this:
    # |           p           |
    # |x|x|x|x|x|x|x|x|...|x|x|
    # |o|x|x|o|x|x|o|x|...|x|o|
    # |x|x|x|x|x|x|x|x|...|x|x|
    # |x|x|x|x|x|x|x|x|...|x|x|
    # |o|x|x|o|x|x|o|x|...|x|o|
    # |x|x|x|x|x|x|x|x|...|x|x|
    # Each successive iteration of the for loop advances the 'o' entries once
    # to the right (for i) and down one row (for j). Once a launch completes,
    # the next wave is enqueued and launched
    print("\nProcessing all remaining zones\n")
    for i in range(0, 3):
        for j in range (0, 3):
            zone_indices = get_active_indices(i, j, stride)
            print("")
            process_indices(zone_indices, args)

    end = time.time()
    print("Time elapsed: %is" % (int(end - start)))

def move_zero_edge_data(save_dir):
    """Moves any data whose containers have an edge with ra < 0 +360 degrees away"""
    global global_tree
    global_tree_filename = save_dir + "/global.json"
    global_tree = QuadTreeNode.from_file(global_tree_filename, leafClass=IDLeaf)

    global height

    print("Processing zones on RA = 0 edge")

    # first move data near RA ~ 0 to RA 360+
    for j in range(1, height - 2):
        # get the left and right edge zones
        zone_id = get_zone_id_from_indices(width-1, j)
        adj_zone_id = get_zone_id_from_indices(0, j)

        # run the fix-overlaps function, but flip the order of the zones
        # to ensure data gets moved from zone_id to adj_zone_id
        fix_overlaps(save_dir, zone_id, [adj_zone_id])

def process_indices(zone_indices, args):

    # set up parallel launches
    pool = Pool(args.jobs)

    # create a partial to store the save directory
    fzo_func = partial(fix_overlaps_from_positions, args.save_dir)

    if args.debug:
        for zone_index in zone_indices:
            #print zone_index
            fzo_func(zone_index)
    else:
        result = pool.imap(fzo_func, zone_indices)
        pool.close()
        pool.join()


if __name__ == "__main__":
    main()
