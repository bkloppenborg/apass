#!/bin/python

# system-level includes
import argparse
import os
from numpy import *
import glob
import numpy as np

# APASS-specific things
import apass
from apass_types import *
from quadtree import *
from quadtree_types import *

# parallel processing
import multiprocessing as mp

# The output looks like this
##  Name    RA(J2000)   raerr  DEC(J2000) decerr nobs  mobs       filt  mag  err
#0020131545  11.198035  0.477 -32.933803  0.396    3   12 12.828  0.931 13.758 13.245 12.593 12.471  0.159  0.164  0.039  0.036  0.088  0.289

def summarize_data(container):

    # A rectangle file should contain data on precisely one star take in
    # multiple photmetric filters.

    # Read in the data. This will be a numpy.narray object, so we can slice
    # and dice it however we would like.

    data = container.data
    dtype={'names': apass.fredbin_col_names,'formats': apass.fredbin_col_types}
    data = np.asarray(container.data, dtype=dtype)

    # return an empty string for containers with no data.
    if len(data) == 0:
        return ""

    # RA and DEC:
    name = data['star'][0]
    ra = average(data['ra'])
    ra_sig = std(data['ra'])
    dec = average(data['dec'])
    dec_sig = std(data['dec'])

    # filter out measurements outside of a 2048 - 100 = 1948 radius from
    # the center of the CCD. (this is for 4096x4096 images only)
    center = 2048
    radius = (2048-100)**2
    radius_2 = radius*radius

    # look up the (x,y) location relative to the center of the CCD
    x = data['ccdx'] - center
    y = data['ccdy'] - center
    ccd_radius_2 = x*x + y*y

    # filter out the measurements outside of that radius
    data = data[ccd_radius_2 < radius_2]

    # If no observations passed the filtering stage, return an empty string
    if len(data) == 0:
        return ""

    # compute the magnitudes
    mags = dict()
    filters = set(data['filter_id'])
    for filter in filters:
        temp = data[data['filter_id'] == filter]
        mag = average(temp['xmag1'])
        mag_sig = std(temp['xmag1'])
        mags[filter] = {'mag': mag, 'sig': mag_sig}

    # compute the number of observations and number of nights that made it through filtering
    num_observations = len(data)
    num_nights = len(set(data['night']))

    # construct the output string
    prefix_fmt = "%010i %10.6f %6.3f %10.6f %6.3f %4i %4i"
    mag_fmt = "%6.3f"
    # Start with the following information
    #  # Name    RA(J2000)   raerr  DEC(J2000) decerr nobs  mobs
    #  # 0020131545  11.198035  0.477 -32.933803  0.396    3   12
    output = (prefix_fmt) % (name, ra, ra_sig, dec, dec_sig, num_nights, num_observations)
    # now follow up with magnitudes and errors. Note that APASS splits them as
    #  (mag1, ..., magN, err1, ..., errN)
    for filter_id in range(0, apass.num_filters):
        mag = 99.999
        if filter_id in mags:
            mag = mags[filter_id]['mag']

        output += " " + (mag_fmt) % (mag)

    for filter_id in range(0, apass.num_filters):
        mag_sig = 99.999
        if filter_id in mags:
            mag_sig = mags[filter_id]['sig']

        output += " " + (mag_fmt) % (mag_sig)

    return output


def zone_to_data(zone_container_filename):
    """Processes all of the rectangles found in zone. Zone should be a valid subdirectory
    of apass.apass_save_dir"""

    zone_id = apass.zone_from_name(zone_container_filename)
    zone_name = apass.name_zone(zone_id)
    print "Processing zone " + zone_name

    # load the zone's tree and data from disk and get the leaves
    zone_json = apass.apass_save_dir + apass.name_zone_json_file(zone_id)
    zone_tree = QuadTreeNode.from_file(zone_json, leafClass=RectLeaf)
    apass.load_zone_data(zone_tree, apass.apass_save_dir)
    leaves = zone_tree.get_leaves()

    # process the data
    outfile_name = apass.apass_save_dir + "/" + zone_name + ".dat"
    with open(outfile_name, 'w')  as outfile:
        for leaf in leaves:
            for container in leaf.containers:
                output = summarize_data(container)

                if len(output) > 0:
                    outfile.write(output + "\n")

def main():

    parser = argparse.ArgumentParser(
        description="Converts a containerized zone into APASS photometric output")
    parser.add_argument('-j','--jobs', type=int, help="Parallel jobs", default=1)
    parser.add_argument('--debug', type=bool, help="Run in debug mode", default=False)
    parser.add_argument('input', nargs='+')
    args = parser.parse_args()

    # run in debug mode
    if args.debug:
        for zonefile in args.input:
            zone_to_data(zonefile)
    # run in production mode
    else:
        pool = mp.Pool(args.jobs)
        result = pool.map_async(zone_to_data, args.input)
        pool.close()
        pool.join()

if __name__ == "__main__":
    main()
