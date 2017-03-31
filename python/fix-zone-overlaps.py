#!/bin/python

# system-level includes
import argparse
import os

# APASS-specific things
from quadtree import *
from quadtree_types import *
import apass
from apass_types import *

def main():

    parser = argparse.ArgumentParser(description='Merges .fred photometric files')

    args = parser.parse_args()

    save_dir = apass.apass_save_dir

    # restore the global tree, we'll need it shortly
    global_tree_file =  save_dir + "/global.json"
    global_tree = QuadTreeNode.from_file(global_tree_file, leafClass=IDLeaf)

    zones = global_tree.get_leaves()

    # iterate over the zones in the global tree. Zones are of type IDLeaf
    for zone in zones:
        # get the zone_id from the node
        zone_id = zone.node_id
        print("Processing Zone: %i" % (zone_id))

        # first get the names of all of the important files
        zone_border_info_file = save_dir + apass.name_zone_border_file(zone_id)
        zone_container_file = save_dir + apass.name_zone_container_file(zone_id)

        # test if the zone has data by checking the presence of the zone's directory
        # if there is no directory, continue
        if not os.path.isfile(zone_container_file):
            print(" no data")
            continue

        # test for the presence of rectangles on the borders of zones by checking for
        # the `border-rects.txt` file. If that doesn't exist, continue
        if not os.path.isfile(zone_border_info_file):
            print(" no zone-border-rects.json file")
            continue

        # get a dictionary containing border rect information
        border_info = apass.load_border_info(zone_border_info_file)

        # The zone exists and there are rectangles on the border. Load the zone's data
        # and the border rectangle data
        zone_file = apass.apass_save_dir + apass.name_zone_json_file(zone_id)
        zone_tree = QuadTreeNode.from_file(zone_file, leafClass=RectLeaf)
        apass.load_zone_data(zone_tree, apass.apass_save_dir)

        # iterate over the border rectangles data
        for key, info in border_info.iteritems():
            # find the RectLeaf node which contains this border rectangle
            print(" checking for overlaps with %s" % (key))
            x,y = info['center']
            node = zone_tree.find_leaf(x,y)
            node_id = node.node_id

            # find the RectContainer which holds this border rectangle
            container = node.get_container(x,y)
            if container is None:
                raise RuntimeError("Could not locate a container holding the border rectangle!")

            # the container data may have been moved, if so move on
            if container.moved:
                continue

            # now that we have a RectContainer, find adjacent containers via the global tree
            # and subsequent zone trees:
            adj_zones = []
            for x,y in container.get_corners():
                # wrap the (x,y) values into 0 <= x <= 360, -90 < y < 90
                x,y = wrap_bounds(x, y, 0, 360, -90, 90)

                # find the adjacent zone, skip it if we've been there recently
                adj_zone = global_tree.find_leaf(x,y)
                if adj_zone == zone or adj_zone in adj_zones:
                    continue
                adj_zones.append(adj_zone)

                # load the adjacent zone tree
                adj_zone_id = adj_zone.node_id
                adj_zone_file = save_dir + apass.name_zone_json_file(adj_zone_id)

                # verify that the adjacent zone contains data
                if not os.path.isfile(adj_zone_file):
                    print("  tried zone %i, but it does not contain data" % (adj_zone_id))
                    continue

                # load the adjacent zone's tree.
                # NOTE: We delay data loading until the last possible moment
                adj_zone_tree = QuadTreeNode.from_file(adj_zone_file, leafClass=RectLeaf)

                # load the adjacent zone's border info file
                adj_border_info = {}
                adj_border_info_filename = apass.apass_save_dir \
                                           + apass.name_zone_border_file(adj_zone_id)
                if os.path.isfile(adj_border_info_filename):
                    adj_border_info = apass.load_border_info(adj_border_info_filename)

                # find the adjacent node
                adj_node = adj_zone_tree.find_leaf(x,y)
                adj_node_id = adj_node.node_id

                # find overlapping containers and merge them. update the
                # border info as necessary
                adj_containers = adj_node.get_overlapping_containers(container, remove=False)

                # if there are no adjacent containers, move on
                if len(adj_containers) == 0:
                    continue

                # load the adjacent zone's data
                apass.load_zone_data(adj_zone_tree, apass.apass_save_dir)

                for adj_container in adj_containers:
                    # generate nice output
                    name = apass.name_rect(adj_container.zone_id, adj_container.node_id,
                                           adj_container.container_id)
                    print("  merging %s " % (name))

                    # load the data and merge.
                    container.merge(adj_container, mark_moved=True)
                    if name in adj_border_info.keys():
                        print("  updating border rect file for zone %i" % (adj_zone_id))
                        del adj_border_info[name]

                # save modifications to the adjacent node and zone
                #print(" updating data in (adjacent) node %s" %
                #      (apass.name_node(adj_zone_id, adj_node_id)))
                apass.save_zone_data(adj_zone_tree, apass.apass_save_dir)
                #print(" updating data in adjacent zone %s" %
                #      (apass.name_zone(adj_zone_id)))
                apass.save_border_info(adj_border_info_filename, adj_border_info)
                QuadTreeNode.to_file(adj_zone_tree, adj_zone_file)

        # save modifications to this node
        #print(" updating data in (primary) node %s" %
        #      (apass.name_node(zone_id, node_id)))
        apass.save_zone_data(zone_tree, apass.apass_save_dir)

        # save modifications to the zone file
        #print(" updating data in (primary) zone %s" %
        #      (apass.name_zone(zone_id)))
        QuadTreeNode.to_file(zone_tree, zone_file)


def wrap_bounds(x, y, x_min, x_max, y_min, y_max):
    # TODO: handle border wrapping
    return x,y

if __name__ == "__main__":
    main()
