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

# numpy array manipulation
import numpy.lib.recfunctions as nprf

# suppress FutureWarning from np.average
import warnings
warnings.simplefilter(action = "ignore", category = FutureWarning)


# The output looks like this
##  Name    RA(J2000)   raerr  DEC(J2000) decerr nobs  mobs       filt  mag  err
#0020131545  11.198035  0.477 -32.933803  0.396    3   12 12.828  0.931 13.758 13.245 12.593 12.471  0.159  0.164  0.039  0.036  0.088  0.289

def summarize_data(container):
    """Parses the measurements contained within a container and generates a
    string summarizing the data."""

    # A rectangle file should contain data on precisely one star take in
    # multiple photmetric filters.

    # Read in the data. This will be a numpy.narray object, so we can slice
    # and dice it however we would like.

    data = container.data
    dtype={'names': apass.fredbin_col_names,'formats': apass.fredbin_col_types}
    data = np.asarray(container.data, dtype=dtype)

    # add a boolean column to serve as a flag for whether or not a given datum should be used
    tmp = np.zeros(len(data))
    data = nprf.append_fields(data, ['use_data'], [tmp], dtypes=[bool])

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

    #
    # apply various filtering to determine when we set the 'use_data' column to true
    filter_ids = set(data['filter_id'])

    # filter by CCD radius these data are always used
    indexes = np.where(ccd_radius_2 < radius_2)
    data['use_data'][indexes] = True

    # if the radius test caused a filter to have fewer than N observations, restore
    # all of the measurements for that filter_id
    for filter_id in filter_ids:
        indexes = np.where(data['filter_id'] == filter_id)
        temp = data[indexes]
        if sum(temp['use_data']) < apass.min_num_observations:
            data['use_data'][indexes] = True

    # If no observations passed the filtering stage, we have a serious problem.
    # Print out a notification and return an empty string. We'll have to let the user
    # figure out what to do.
    if sum(data['use_data']) == 0:
        zone_id = data['zone_id'][0]
        node_id = data['node_id'][0]
        container_id = data['container_id'][0]
        print("WARNING: No measurements passed filtering for zone %i node %i container %i" % (zone_id, node_id, container_id))
        return ""

    # compute the magnitudes
    mags = dict()
    for filter_id in filter_ids:
        # extract known-good measurements for this filter
        temp = data[(data['filter_id'] == filter_id) & (data['use_data'] == True)]
        temp = data[indexes]
        num_obs = len(temp)
        # compute the average and standard deviation for the magnitude.
        # NOTE: Using standard error propigation methods, e.g.
        #        mag_sig = sqrt(sum(error_i^2) / N
        #       double-counts the nightly photometric (Poisson) error.
        #       Arne indicates we should instead use std(mag) to avoid this problem
        mag = average(temp['xmag1'])
        mag_sig = std(temp['xmag1'])
        mags[filter_id] = {'mag': mag, 'sig': mag_sig}

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

    out_mags = []
    out_mag_sigs = []
    for filter_id in range(1, apass.num_filters + 1):

        # The existing APASS pipeline computes the six output filters
        # as follows:
        # index | out <- in
        #   1 <- 3
        #   2 <- (2 - 3)
        #   3 <- 2
        #   4 <- 8
        #   5 <- 9
        #   6 <- 10
        # we mirror this here
        if filter_id == 1:
            mag, mag_sig = read_mags(mags, 3)
        elif filter_id == 2:
            mag2, mag2_sig = read_mags(mags, 2)
            mag3, mag3_sig = read_mags(mags, 3)

            if mag2 == 99.999 or mag3 == 99.999:
                mag = 99.999
                mag_sig = 99.999
            else:
                mag = mag2 - mag3
                mag_sig = sqrt(mag2_sig**2 + mag3_sig**2)
        elif filter_id == 3:
            mag, mag_sig = read_mags(mags, 2)
        elif filter_id == 4:
            mag, mag_sig = read_mags(mags, 8)
        elif filter_id == 5:
            mag, mag_sig = read_mags(mags, 9)
        elif filter_id == 6:
            mag, mag_sig = read_mags(mags, 10)

        out_mags.append(mag)
        out_mag_sigs.append(mag_sig)

    # write out the magnitudes
    for filter_id in range(0, apass.num_filters):
        output += " " + (mag_fmt) % (out_mags[filter_id])

    # and the corresponding uncertainties
    for filter_id in range(0, apass.num_filters):
        output += " " + (mag_fmt) % (out_mag_sigs[filter_id])

    return output

def read_mags(mags_dict, filter_id):
    """Reads the magnitude and error from the dictionary for the specified filter."""
    mag = 99.999
    mag_sig = 99.999

    if filter_id in mags_dict:
        mag = mags_dict[filter_id]['mag']
        mag_sig = mags_dict[filter_id]['sig']

    return [mag, mag_sig]

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
    parser.add_argument('-j','--jobs', type=int, help="Parallel jobs", default=4)
    parser.add_argument('input', nargs='+')
    parser.add_argument('--debug', default=False, action='store_true',
                        help="Run in debug mode")
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
