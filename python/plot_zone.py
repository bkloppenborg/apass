#!/bin/python

import argparse
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.patches as patches

import apass
from apass_types import *
from quadtree import *
from quadtree_types import *

# file I/O
import fred

def plot_containers(leaf, axes):
    """Adds patches corresponding to container borders"""

    for container in leaf.containers:
        rect = container.rect
        x = rect.x_min
        y = rect.y_min
        dx = rect.x_max - x
        dy = rect.y_max - y

        axes.add_patch(patches.Rectangle((x,y), dx, dy, color="red", fill=False))

def main():
    """Plots all data in an APASS zone and optionally displays container boundaries."""

    parser = argparse.ArgumentParser(description='Plots a zone file and its data')
    parser.add_argument('input', nargs='+', help="Zone JSON files to be plotted")
    parser.add_argument('--show-containers', default=False, action="store_true",
                        help="Plot borders for containers")

    args = parser.parse_args()
    save_dir = os.path.dirname(os.path.realpath(args.input[0])) + "/"

    inputs = args.input

    for filename in inputs:
        zone_id = apass.zone_from_name(filename)

        zone_name = apass.name_zone(zone_id)
        print "Plotting zone " + zone_name

        # load the original zone data. Note, we don't restore it to the tree
        zone_data_file = save_dir + apass.name_zone_file(zone_id)
        zone_data = fred.read_fredbin(zone_data_file)
        print("Zone file has " + str(zone_data.size) + " entries")

        # load the containerized zone data
        zone_container_file = save_dir + apass.name_zone_container_file(zone_id)
        zone_container_data = fred.read_fredbin(zone_container_file)
        print("Zone container file has " + str(zone_container_data.size) + " entries")

        # load the zone's tree
        zone_json = save_dir + apass.name_zone_json_file(zone_id)
        zone_tree = QuadTreeNode.from_file(zone_json, leafClass=RectLeaf)
        leaves = zone_tree.get_leaves()

        total_containers = 0
        for leaf in leaves:
            total_containers += len(leaf.containers)

        print("Zone contains a total of " + str(total_containers) + " containers")

        fig, axes = plt.subplots(1)

        # plot the zone container data
        for leaf in leaves:
            rect = leaf.rect
            x = rect.x_min
            y = rect.y_min
            dx = rect.x_max - x
            dy = rect.y_max - y

            axes.add_patch(patches.Rectangle((x,y), dx, dy, fill=False))

            if args.show_containers:
                plot_containers(leaf, axes)

        # plot the data
        dec = zone_data['dec']
        ra  = zone_data['ra']
        plt.scatter(ra, dec)

        dec = zone_container_data['dec']
        ra  = zone_container_data['ra']
        plt.scatter(ra, dec, color="green")

        plt.show()


if __name__ == "__main__":
    main()
