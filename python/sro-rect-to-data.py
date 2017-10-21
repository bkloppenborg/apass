#!/bin/python

# system-level includes
import argparse
import os
from numpy import *
import glob
import numpy as np
import time
import itertools
import networkx as nx

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

# the maximum radius for which the field flattener works correctly
sro_ccd_x_center = 2048 # pixels
sro_ccd_y_center = 2048 # pixels
sro_max_ccd_radius = 3072.0 / 2
sro_filter_ids = [2, 3, 8, 9, 10] # B, V, sg, sr, si
sro_num_filters = len(sro_filter_ids) + 1 # we write out V, (B-V), B, sg, sr, si)
sro_min_num_observations = 3

# The output looks like this
##  Field    RA(J2000)   raerr  DEC(J2000) decerr nobs  mobs       filt  mag  err
#0020131545  11.198035  0.477 -32.933803  0.396    3   12 12.828  0.931 13.758 13.245 12.593 12.471  0.159  0.164  0.039  0.036  0.088  0.289

def filter_by_ccd_radius(data, x_center, y_center, max_ccd_radius):
    """Filters """

    # filter out measurements outside of a specific radius from the center
    # of the CCD
    radius = (max_ccd_radius)**2
    radius_2 = radius*radius

    # Compute the (x,y) location of the stars relative to the center of the CCD
    x = data['ccdx'] - x_center
    y = data['ccdy'] - y_center
    ccd_radius_2 = x*x + y*y

    # filter by CCD radius these data are always used
    indexes = np.where(ccd_radius_2 < radius_2)
    data['use_data'][indexes] = True

    return data

def average_by_field(container):
    """Parses the measurements contained within a container and averages the photometry
    for each field number.

    This function returns:
    1. A list of string containing the averaged data
    2. A list of field IDs stored in this container.
    """

    # A container should store data on precisely one star, but that data will
    # be taken through multiple photometric filters and potentially originate from
    # multiple fields (telescope pointings). In this function, we average
    # the photometry for each field and generate one output line for each field.

    # Read in the data. This will be a numpy.narray object, so we can slice
    # and dice it however we would like.
    data = container.data
    dtype={'names': apass.fredbin_col_names,'formats': apass.fredbin_col_types}
    data = np.asarray(container.data, dtype=dtype)

    # add a boolean column to serve as a flag for whether or not a given datum should be used
    tmp = np.zeros(len(data))
    data = nprf.append_fields(data, ['use_data'], [tmp], dtypes=[bool])

    # if the container holds no data, immediately return an empty string
    if len(data) == 0 or container.moved == True:
        return ""

    # Compute the average RA and DEC using data from all measurements
    ra = average(data['ra'])
    ra_sig = std(data['ra'])
    dec = average(data['dec'])
    dec_sig = std(data['dec'])

    ##
    # Filtering Stages
    ##
    data = filter_by_ccd_radius(data, sro_ccd_x_center, sro_ccd_y_center, sro_max_ccd_radius)

    field_ids  = set(data['field'])
    filter_ids = set(data['filter_id'])

    # if the radius test caused a filter to have fewer than N observations, restore
    # all of the measurements for that filter_id
    for filter_id in filter_ids:
        indexes = np.where(data['filter_id'] == filter_id)
        temp = data[indexes]
        if sum(temp['use_data']) < sro_min_num_observations:
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

    # Now average the photometry on a per-field basis:
    output = []
    for field_id in field_ids:
        indexes = np.where((data['field'] == field_id))
        t_data = data[indexes]

        mags = average_magnitudes(t_data, filter_ids)

        # compute the number of observations and number of nights that made it through filtering
        num_observations = len(t_data)
        num_nights = len(set(t_data['night']))

        # construct the output string
        prefix_fmt = "%25s %10.6f %6.3f %10.6f %6.3f %4i %4i"
        mag_fmt = "%6.3f"
        # Start with the following information
        #  # Field    RA(J2000)   raerr  DEC(J2000) decerr nobs  mobs
        #  # 0020131545  11.198035  0.477 -32.933803  0.396    3   12
        line = (prefix_fmt) % (field_id, ra, ra_sig, dec, dec_sig, num_nights, num_observations)
        # now follow up with magnitudes and errors. Note that APASS splits them as
        #  (mag1, ..., magN, err1, ..., errN)

        # Loop over the filters and apply any custom reorganization to them.
        out_mags = []
        out_mag_sigs = []
        for filter_id in range(1, sro_num_filters + 1):

            # SRO filters are computed as follows
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
        for filter_id in range(0, sro_num_filters):
            line += " " + (mag_fmt) % (out_mags[filter_id])

        # and the corresponding uncertainties
        for filter_id in range(0, sro_num_filters):
            line += " " + (mag_fmt) % (out_mag_sigs[filter_id])

        output.append(line)

    # All done, return the string.
    if len(output) != len(field_ids):
        raise "Field averages and number of fields do not match!"

    return output, field_ids

def average_magnitudes(data, filter_ids):

    # compute the magnitudes
    mags = dict()
    for filter_id in filter_ids:
        # extract known-good measurements for this filter
        indexes = np.where((data['filter_id'] == filter_id) & (data['use_data'] == True))
        temp = data[indexes]
        num_obs = len(temp)

        if num_obs > 0:

            # compute the average and standard deviation for the magnitude.
            # NOTE: Using standard error propigation methods, e.g.
            #        mag_sig = sqrt(sum(error_i^2) / N
            #       double-counts the nightly photometric (Poisson) error.
            #       Arne indicates we should instead use std(mag) to avoid this problem
            # If there is only one observation, copy the error.
            mag = average(temp['xmag1'])
            if num_obs > 1:
                mag_sig = std(temp['xmag1'])
            else:
                mag_sig = temp['xerr1']

        else:
            mag = 99.999
            mag_sig = 99.999

        mags[filter_id] = {'mag': mag, 'sig': mag_sig, 'num': num_obs}

    return mags

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

    # create a graph data structure for this zone
    G = nx.Graph()

    # load the zone's tree and data from disk and get the leaves
    zone_json = apass.apass_save_dir + apass.name_zone_json_file(zone_id)
    zone_tree = QuadTreeNode.from_file(zone_json, leafClass=RectLeaf)
    apass.load_zone_data(zone_tree, apass.apass_save_dir)
    leaves = zone_tree.get_leaves()

    # process the data
    line_no = 0
    dat_filename = apass.apass_save_dir + "/" + zone_name + ".dat"
    with open(dat_filename, 'w') as dat_file:
        for leaf in leaves:
            for container in leaf.containers:
                output_lines, field_ids = average_by_field(container)

                # save the text data to file, noting what lines it gets
                # written to
                line_numbers = []
                for line in output_lines:
                    dat_file.write(line + "\n")
                    line_numbers.append(line_no)
                    line_no += 1

                # for stars that have more than one field, we add this information
                # to the graph data structure to recover the data in the next stage
                # of the pipeline
                if len(field_ids) > 1:
                    # generate all possible combinations of the line numbers and field IDs
                    # these express the edges in the graph
                    line_pairs = list(itertools.combinations(line_numbers, 2))
                    field_pairs = list(itertools.combinations(field_ids, 2))

                    for line_pair, field_pair in zip(line_pairs, field_pairs):
                        src_line = line_pair[0]
                        dst_line = line_pair[1]
                        src_field = field_pair[0]
                        dst_field = field_pair[1]

                        try:
                            edge = G[src_field][dst_field]
                        except:
                            G.add_edge(src_field, dst_field, line_ids=[], weight=0)
                            edge = G[src_field][dst_field]

                        edge['line_ids'].append((src_line, dst_line))
                        edge['weight'] += 1

    # save the pickle file
    graph_filename = apass.apass_save_dir + "/" + zone_name + ".p"
    nx.write_gpickle(G, graph_filename)


def main():

    parser = argparse.ArgumentParser(
        description="Converts a containerized zone into APASS photometric output")
    parser.add_argument('-j','--jobs', type=int, help="Parallel jobs", default=4)
    parser.add_argument('input', nargs='+')
    parser.add_argument('--debug', default=False, action='store_true',
                        help="Run in debug mode")
    args = parser.parse_args()

    start = time.time()
    # run in debug mode
    if args.debug:
        for zonefile in args.input:
            zone_to_data(zonefile)
    # run in production mode
    else:
        pool = mp.Pool(args.jobs)
        result = pool.imap(zone_to_data, args.input)
        pool.close()
        pool.join()

    end = time.time()
    print("Time elapsed: %is" % (int(end - start)))

if __name__ == "__main__":
    main()
