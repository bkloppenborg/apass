#!/usr/bin/python
# This script finds zones that have not been completely processed

# system includes
import argparse
import glob
import os

# local includes
from quadtree import *
from quadtree_types import *
from apass import name_zone, north_zone_id, south_zone_id

expected_extensions = [
    ".fredbin",
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
    poles_checked = [False, False]
    leaves = tree.get_leaves()
    for leaf in leaves:
        zone_id = leaf.node_id
        zone_name = name_zone(zone_id)

        # The north/south zone IDs appear multiple times in the tree.
        # Ensure we only visit these once.
        if zone_id == south_zone_id or zone_id == north_zone_id:
            if poles_checked[zone_id]:
                continue
            else:
                poles_checked[zone_id] = True

        # determine which, if any, files are missing for this zone:
        missing_extensions = []
        for extension in expected_extensions:
            filename = save_dir + zone_name + extension

            if filename in saved_files:
                saved_files.remove(filename)
            else:
                missing_extensions.append(extension)

        if len(missing_extensions) > 0:
            print str(zone_id) + ": " + ", ".join(missing_extensions)

if __name__ == "__main__":
    main()
