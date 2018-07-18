#!/usr/bin/python

# system includes
import argparse
import os
import numpy as np

# locals
from quadtree import *
from quadtree_types import *
from apass_types import IDLeaf, RectLeaf
import apass
import zone
import fred

# import types from rect_to_dat
from rect_to_dat import filter_config_data # class
# import functions from rect_to_dat
from rect_to_dat import apply_filters, compute_weights
from rect_to_dat import apass_ccd_x_center, apass_ccd_y_center, apass_max_ccd_radius
from rect_to_dat import apass_min_num_observations

def data_from_coordinate_pair(save_dir, x,y):

    # load the global tree, find the zone_id to which this star belongs
    global_tree_file = save_dir + "/global.json"
    global_tree      = QuadTreeNode.from_file(global_tree_file, leafClass = IDLeaf)
    leaf             = global_tree.find_leaf(x,y)
    zone_id          = leaf.node_id

    # load the zone, extract the star's data
    zone_name = apass.name_zone(zone_id)

    # load the zone's tree and data from disk and get the leaves
    zone_json = save_dir + apass.name_zone_json_file(zone_id)
    zone_tree = QuadTreeNode.from_file(zone_json, leafClass=RectLeaf)
    zone.load_zone_data(zone_tree, save_dir)

    # find the leaf containing the point of interest
    leaf = zone_tree.find_leaf(x,y)

    # export the data
    data = None
    for container in leaf.containers:
        if container.contains(x,y):
            data = fred.to_fredbin(container.data)

    return data

def data_from_unique_id(save_dir, zone_id, node_id, container_id):

    zone_filename = save_dir + apass.name_zone(zone_id) + "-container.fredbin"
    zone_data = fred.read_fredbin(zone_filename)

    indices = np.where((zone_data['node_id'] == node_id) &
                       (zone_data['container_id'] == container_id))

    data = zone_data[indices]

    return data

def main():

    # output globals
    parser = argparse.ArgumentParser(description='Dumps all data for a given star.')
    parser.add_argument('save_dir', help="Directory where save files can be found")
    parser.add_argument('coords', nargs='+',
                        help="Coordinates in either a (x,y) pair, " + \
                        "or (zone_id, node_id, container_id) triplet")

    # parse the command line arguments
    args = parser.parse_args()

    coords = args.coords
    num_coords = len(coords)
    if num_coords == 2:
        coords = map(float, args.coords)
    elif num_coords == 3:
        coords = map(int, args.coords)
    else:
        print("Coordinates must be a pair or triplet. See -h for information.")
        quit()

    save_dir = os.path.abspath(args.save_dir) + "/"

    # configure the filters
    filter_config = filter_config_data()
    filter_config.load_apass_defaults()

    # export the data
    zone_id, node_id, container_id = [None, None, None]
    data = None
    if num_coords == 2:
        x,y = coords
        data = data_from_coordinate_pair(save_dir, x, y)
        if data is None:
            print("Could not locate star.")
            quit()
        zone_id      = data['zone_id'][0]
        node_id      = data['node_id'][0]
        container_id = data['container_id'][0]
    elif num_coords == 3:
        zone_id, node_id, container_id = coords
        data = data_from_unique_id(save_dir, zone_id, node_id, container_id)
        if data is None:
            print("Could not locate star.")
            quit()

    # convert the fredbin to a freddat, apply filters, and compute weights
    # as would be performed in the rect_to_dat code.
    data = fred.to_freddat(data)
    data = apply_filters(data, filter_config)
    data = compute_weights(data)

    data.sort(order='filter_id')

    # write the star's data to disk
    out_filename = save_dir + 'star-' + str(zone_id) + '-' + \
                   str(node_id) + '-' + str(container_id)
    fred.write_freddat_txt(data, out_filename)


if __name__ == "__main__":
    main()
