#!/bin/python

# system includes
import argparse
import os.path

# local includes
from quadtree import *
from quadtree_types import *

def main():
    """Exports zone metadata to a CSV file"""

    parser = argparse.ArgumentParser(description='Exports zone information to CSV file.')
    parser.add_argument('tree_file', \
                        help="The global.json zone definition file.")

    args = parser.parse_args()
    tree_file = args.tree_file
    tree_file_csv = os.path.splitext(tree_file)[0] + ".csv"

    # restore the tree. make-zones.py writes out leaves of type IDLeaf
    tree = QuadTreeNode.from_file(tree_file, leafClass=IDLeaf)

    # get the leaves and build a list containing zone_id, ra, dec
    leaves = tree.get_leaves()

    seen_zone_id0 = False
    seen_zone_id1 = False
    zone_infos = []
    for leaf in leaves:
        zone_id = leaf.node_id
        x_min = leaf.rect.x_min
        x_max = leaf.rect.x_max
        y_min = leaf.rect.y_min
        y_max = leaf.rect.y_max

        # be sure to write out only one instance of zone_id 0 and 1
        if zone_id == 0:
            if seen_zone_id0:
                continue
            else:
                seen_zone_id0 = True

        if zone_id == 1:
            if seen_zone_id1:
                continue
            else:
                seen_zone_id1 = True

        zone_infos.append([zone_id, x_min, x_max, y_min, y_max])

    # sort by zone_id
    zone_infos.sort(key=lambda x: x[0])

    # remove duplicate zone_id = 0, 1 entries

    with open(tree_file_csv, 'w') as outfile:
        outfile.write("zone_id, ra_min, ra_max, dec_min, dec_max \n")
        for zone_info in zone_infos:
            zone_info = map(str, zone_info)
            outfile.write(','.join(zone_info) + "\n")


if __name__ == "__main__":
    main()
