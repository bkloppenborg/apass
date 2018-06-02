#!/usr/bin/python

# system includes
import argparse
import os

# locals
from quadtree import *
from quadtree_types import *
from apass_types import IDLeaf
from apass import name_zone

def main():

    parser = argparse.ArgumentParser(description='Identifies the state of the pipeline for each zone')
    parser.add_argument('save_dir', help="Directory where save files can be found")
    parser.add_argument('coords', nargs='+',
                        help="Coordinates in either a (x,y) pair, or (x_min, y_min, x_max, y_max) quad")

    # parse the command line arguments
    args = parser.parse_args()

    coords = map(float, args.coords)
    num_coords = len(coords)
    if not (num_coords == 2 or num_coords == 4):
        print("Coordinates must either be a pair or a quad. See -h")
        quit()

    save_dir = os.path.abspath(args.save_dir) + "/"
    tree_file   = save_dir + "/global.json"

    # load the zone quadtree
    tree = QuadTreeNode.from_file(tree_file, leafClass=IDLeaf)

    # execute the search
    zone_ids = []
    if num_coords == 2:
        # simple (x,y) pair, simply find the leaf
        x,y = coords
        leaf = tree.find_leaf(x,y)
        zone_ids.append(leaf.node_id)
    elif num_coords == 4:
        print("Not implemented")

    for zone_id in zone_ids:
        print(name_zone(zone_id))

if __name__ == "__main__":
    main()
