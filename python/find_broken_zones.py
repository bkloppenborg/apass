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

# these are file extensions we expect to find in the save directory
expected_extensions = [
    ".fredbin",              # stage 1
    "-contrib.txt",          # stage 1
    "-container.fredbin",    # stage 2
    "-border-rects.json",    # stage 2
    ".dat"                   # stage 3
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

    # init a structure to store stage -> zone mapping information
    stages = dict()
    for extension in expected_extensions:
        stages[extension] = []

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
                stages[extension].append(filename)

        if len(missing_extensions) > 0:
            # print a message on the console
            print str(zone_id) + ": " + ", ".join(missing_extensions)


    missing_file = save_dir + '/missing.txt'
    for extension in expected_extensions:
        missing_files = stages[extension]
        if len(missing_files) > 0:
            print("The earliest files missing are %s." % (extension))
            print("Modify the 'missing.txt' file for the previous stage input " +
                  "files and run the corresponding stage of the pipeline")
            with open(missing_file, 'w') as outfile:
                for filename in missing_files:
                    outfile.write(filename + "\n")
            break

if __name__ == "__main__":
    main()
