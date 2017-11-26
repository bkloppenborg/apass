#!/usr/bin/python

# system includes
import argparse
import sys, os
from multiprocessing import Pool

# for plots
import numpy as np
import matplotlib.pyplot as plt

# local includes
import apass
from apass_types import *
from quadtree import QuadTreeNode
from quadtree_types import *

# custom modules
sys.path.append(os.path.join(sys.path[0],'modules', 'FileLock', 'filelock'))
from filelock import FileLock

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


def get_adjacent_zones(i, j):
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

    # polar zones, return an entire row
    if j == 0 or j == height - 1:
        for k in range(0, width):
            if k == i:
                continue

            adj_zone_indices.append((k, j))

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

    adj_zones = []
    for k,l in adj_zone_indices:
        adj_zone = get_zone_from_indices(k, l)
        adj_zones.append(adj_zone)
    return adj_zones

def get_zone_from_indices(i,j):
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

    return get_zone_from_coordinates(x,y)

def get_zone_from_coordinates(x,y):
    """Returns a reference to a zone in the global tree given the specified
    (x,y) coordinates.

    Requires the global_tree is loaded into the global namespace."""
    global global_tree
    return global_tree.find_leaf(x,y)

def load_zone(zone_id):
    """Loads the tree, data, and border info file for the specified zone.
    Returns this data as a dictionary keyed as follows:
     * json_filename
     * container_filename
     * border_info_filename
     * border_info
     * tree
     * loaded
     * lock (a FileLock instance)"""

    #print(" Loading zone %i" % (zone_id))

    zone_dict = dict()
    zone_tree = None
    zone_border_info = None
    load_succeeded = True
    lock = None

    zone_json_file        = apass.apass_save_dir + apass.name_zone_json_file(zone_id)
    zone_container_file   = apass.apass_save_dir + apass.name_zone_container_file(zone_id)
    zone_border_info_file = apass.apass_save_dir + apass.name_zone_border_file(zone_id)

    lock = FileLock(zone_container_file)
    lock.acquire()

    # verify the necessary files exist, if not bail
    if not os.path.isfile(zone_json_file):
#        print(" WARNING: JSON file missing for zone %i" % (zone_id))
        load_succeeded = False
    if not os.path.isfile(zone_border_info_file):
#        print(" WARNING: Border Info file missing for zone %i" % (zone_id))
        load_succeeded = False
    if not os.path.isfile(zone_container_file):
#        print(" WARNING: Container data file missing for zone %i" % (zone_id))
        load_succeeded = False

    if load_succeeded:
        zone_border_info = apass.load_border_info(zone_border_info_file)
        zone_tree = QuadTreeNode.from_file(zone_json_file, leafClass=RectLeaf)
        apass.load_zone_data(zone_tree, apass.apass_save_dir)

    zone_dict['json_filename']        = zone_json_file
    zone_dict['container_filename']   = zone_container_file
    zone_dict['border_info_filename'] = zone_border_info_file
    zone_dict['tree']                 = zone_tree
    zone_dict['border_info']          = zone_border_info
    zone_dict['loaded']               = load_succeeded
    zone_dict['lock']                 = lock
    return zone_dict

def save_zone(zone_dict):
    """Saves the zone data to disk. The zone_dict format matches the
    format found in load_zone above."""

    json_file        = zone_dict['json_filename']
    container_file   = zone_dict['container_filename']
    border_info_file = zone_dict['border_info_filename']
    tree             = zone_dict['tree']
    border_info      = zone_dict['border_info']
    load_succeeded   = zone_dict['loaded']
    lock             = zone_dict['lock']

    zone_id = apass.zone_from_name(zone_container_file)
    #print("Saving zone %i" % (zone_id))

    # save the data to disk
    apass.save_zone_data(tree, apass.apass_save_dir)
    apass.save_border_info(border_info_file, border_info)
    QuadTreeNode.to_file(tree, json_file)

    # release the lock
    lock.release()

def fix_zone_overlaps(zone_index):
    """Fixes zone overlaps between adjacent nodes.

    Requires global_tree_filename is loaded into the global namespace."""

    # load the global tree
    global global_tree_filename
    global global_tree

    # load the global tree
    global_tree = QuadTreeNode.from_file(global_tree_filename, leafClass=IDLeaf)

    # Get this zone and its ID
    i,j = zone_index
    zone_id = get_zone_from_indices(i,j).node_id
    print("Processing Zone: %i" % (zone_id))

    # Load the zone
    zone_dict = load_zone(zone_id)
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
    for adj_zone in get_adjacent_zones(i, j):
        z_id = adj_zone.node_id
        adj_zone_dicts[z_id] = load_zone(z_id)

    # iterate over this zone's border rectangles and merge overlapping entries
    for key, info in border_infos.iteritems():
        #print(" checking for overlaps with %s " % (key))

        # Find this zone's node and container which store this border rectangle
        x,y = info['center']
        node = tree.find_leaf(x,y)
        node_id = node.node_id
        container = node.get_container(x,y)

        # Get the coordinates for the corners of this container
        corners = container.get_corners()

        # find adjacent nodes which overlap with this container
        for c_x, c_y in corners:
            c_x, c_y = wrap_bounds(c_x, c_y)

            adj_zone = global_tree.find_leaf(c_x, c_y)
            adj_zone_id = adj_zone.node_id

            # skip corners that are within this zone
            if adj_zone_id == zone_id:
                continue

            # skip zones that didn't load
            if adj_zone_dicts[adj_zone_id]['loaded'] == False:
                continue

            # pull out some variables for easier referencing
            adj_tree = adj_zone_dicts[adj_zone_id]['tree']
            adj_border_infos = adj_zone_dicts[adj_zone_id]['border_info']

            # find the adjacent node and any adjacent containers
            adj_node = adj_tree.find_leaf(c_x, c_y)
            adj_containers = adj_node.get_overlapping_containers(container, remove=False)

            # move any adjacent containers into this container and mark the moved
            for adj_container in adj_containers:
                # generate nice output
                name = apass.name_container(adj_container.zone_id,
                                            adj_container.node_id,
                                            adj_container.container_id)
                print(" merging %s " % (name))
                container.merge(adj_container, mark_moved=True)
                if name in adj_border_infos.keys():
                    #print(" updating border rect file for zone %i" % (adj_zone_id))
                    del adj_border_infos[name]

    # All done with the merge. Write out data for all files.
    save_zone(zone_dict)
    for key, value in adj_zone_dicts.iteritems():
        save_zone(value)

def wrap_bounds(ra, dec):
    """Wraps (ra,dec) coordinates at boundaries."""
    if dec < -90 or dec > 90:
        ra = (ra + 180) % 360
        dec = 90 - (dec + 90) % 180
    elif ra < 0 or ra > 360:
        ra = (ra + 360) % 360

    return ra, dec


def main():
    parser = argparse.ArgumentParser(
        description='Resolves instances where star data spans a zone boundary')
    parser.add_argument('-j','--jobs', type=int, help="Parallel jobs", default=4)
    parser.add_argument('--debug', default=False, action='store_true',
                        help="Run in debug mode")
    parser.set_defaults(jobs=1)

    args = parser.parse_args()

    global global_tree_filename
    global_tree_filename = apass.apass_save_dir + "/global.json"

    # determine the size of the image
    global width
    global height
    width = 2**(apass.global_depth)
    height = width

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


def process_indices(zone_indices, args):

    # set up parallel launches
    pool = Pool(args.jobs)

    if args.debug:
        for zone_index in zone_indices:
            #print zone_index
            fix_zone_overlaps(zone_index)
    else:
        result = pool.imap(fix_zone_overlaps, zone_indices)
        pool.close()
        pool.join()


if __name__ == "__main__":
    main()