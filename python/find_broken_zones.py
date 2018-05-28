#!/usr/bin/python
# This script finds zones that have not been completely processed

# system includes
import argparse
import glob
import os

# local includes
from quadtree import *
from quadtree_types import *
from apass import name_zone

expected_extensions = [
    ".fredbin"
    "-contrib.txt",
    "-container.fredbin",
    "-border-rects.json",
    ".dat"
]


def main():

    parser = argparse.ArgumentParser(description='Parses .fred files into zone .fredbin files')
    parser.add_argument('save_dir', help="Directory where save files can be found")

    # parse the command line arguments and start timing the script
    args = parser.parse_args()
    save_dir = os.path.abspath(args.save_dir) + "/"

    tree_file   = save_dir + "/global.json"
    output_file = save_dir + "/broken_zones.log"

    # get a list of all of the files in the save directory
    saved_files = glob.glob(save_dir + "z*")

    # load the zone quadtree
    tree = QuadTreeNode.from_file(tree_file, leafClass=IDLeaf)

    # find the broken zones
    broken_zones = []
    zone_0_checked = False
    zone_1_checked = False
    leaves = tree.get_leaves()
    for leaf in leaves:
        zone_id = leaf.node_id

        if zone_id == 0 and zone_0_checked:
            continue
        else:
            zone_0_checked = True

        if zone_id == 1 and zone_1_checked:
            continue
        else:
            zone_1_checked = True

        zone_name = name_zone(zone_id)

        for extension in expected_extensions:
            filename = save_dir + zone_name + extension

            if filename in saved_files:
                saved_files.remove(filename)
            else:
                broken_zones.append(zone_name)
                print zone_name
                break




if __name__ == "__main__":
    main()
