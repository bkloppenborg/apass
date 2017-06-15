#!/bin/python

import argparse
import numpy as np
from multiprocessing import Pool

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import apass
from apass_types import *
from quadtree import *
from quadtree_types import *
import glob
import sys

# temporary for testing
import random

def coords_to_pixel(ra, dec, width, height):
    x = int(ra / 360 * width)
    y = int(-dec / 180 * height) + height / 2

    return [x,y]

def zone_measurement_stats(zone_file):

    print("Processing zone file %s" % (zone_file) )
    zone_tree = QuadTreeNode.from_file(zone_file, leafClass=RectLeaf)

    # compute the (x,y) location of the top-level node in the plot
    rect = zone_tree.rect
    ra = (rect.x_max + rect.x_min) / 2
    dec = (rect.y_max + rect.y_min) / 2

    # determine the min/ave/max measurements.
    num_containers = 0
    min_measurements = sys.maxint
    ave_measurements = 0
    max_measurements = 0
    leaves = zone_tree.get_leaves()
    for leaf in leaves:
        for container in leaf.containers:
            num_data = container.num_data

            if(num_data < min_measurements):
                min_measurements = num_data
            elif(num_data > max_measurements):
                max_measurements = num_data

            ave_measurements += num_data
            num_containers += 1

    # compute the average
    if(num_containers > 0):
        ave_measurements = ave_measurements / num_containers

    return [ra, dec, min_measurements, ave_measurements, max_measurements]


def main():
    """Creates a coverage map which shows the number of observations in each
    node in each zone."""

    # the coverage map is (internally) represented as a (square) Numpy NDArray
    # whose dimensions are defined by the number of times the sphere is sudivided.

    parser = argparse.ArgumentParser(description='Produces a coverage map with data from all zones')
    parser.add_argument('-j', '--jobs', type=int, help="Parallel Jobs", default=4)
    parser.add_argument('--debug', default=False, action='store_true',
                        help="Run in debug mode")
    args = parser.parse_args()

    num_splits = apass.global_depth

    width = 2**num_splits
    height = width
    print("Image size: %i %i" % (width, height))

    # compute statistics in parallel
    zone_files = glob.glob(apass.apass_save_dir + "/*-zone.json")
    results = []
    if args.debug:
        for filename in zone_files:
            results.append(zone_measurement_stats(filename))
    else:
        # generate a pool of threads to process the input
        pool = Pool(args.jobs)

        # farm out the work
        results = pool.imap(zone_measurement_stats, zone_files)
        pool.close()
        pool.join()


    # we will generate three images showing the min, average, and maximum
    # number of measurements.
    min_image = np.zeros((width, height), dtype=float)
    ave_image = np.zeros((width, height), dtype=float)
    max_image = np.zeros((width, height), dtype=float)

    # Iterate over the zone files computing the min, ave, and max number of
    # measurements
    for result in results:
        # unpack the result
        [ra, dec, min_measurements, ave_measurements, max_measurements] = result

        [x,y] = coords_to_pixel(ra, dec, width, height)
        #print("center: %3.4f %3.4f -> %i %i" % (ra, dec, x, y))

        # update the image
        min_image[y,x] = min_measurements
        ave_image[y,x] = ave_measurements
        max_image[y,x] = max_measurements


    vmin = min(min_image.min(), max_image.min())
    vmax = max(min_image.max(), max_image.max())

    fig = plt.figure(figsize=(9,9))
    gs = gridspec.GridSpec(3, 1)
    min_ax = plt.subplot(gs[0, 0])
    im = min_ax.imshow(min_image, extent=[0, 360, -90, 90])
    min_ax.set_title("Minimum number of measurements")
    plt.colorbar(im)

    ave_ax = plt.subplot(gs[1, 0], sharex=min_ax, sharey=min_ax)
    im = ave_ax.imshow(ave_image, extent=[0, 360, -90, 90], vmin=0, vmax=50)
    ave_ax.set_title("Average number of measurements")
    plt.colorbar(im)

    max_ax = plt.subplot(gs[2, 0], sharex=min_ax, sharey=min_ax)
    im = max_ax.imshow(max_image, extent=[0, 360, -90, 90], vmin=0, vmax=100)
    max_ax.set_title("Maximum number of measurements")
    plt.colorbar(im)

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
