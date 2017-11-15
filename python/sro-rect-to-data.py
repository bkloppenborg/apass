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

def filter_by_ccd_radius(data, x_center, y_center, max_ccd_radius,
                         min_num_observations = 3):
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

    # if there are fewer than 'min_num_observations' in a given photometric filter,
    # restore the observations.
    filter_ids = set(data['filter_id'])
    for filter_id in filter_ids:
        indexes = np.where(data['filter_id'] == filter_id)
        temp = data[indexes]
        if sum(temp['use_data']) < min_num_observations:
            data['use_data'][indexes] = True

    return data

def average_by_field(container):
    """Parses the measurements within a container and averages the data by field.

    This function returns a nested dictionary whose keys are the field IDs and values
    correspond to the output of average_container below.
    """

    output = []

    # Get a list of the unique field IDs in the container
    dtype={'names': apass.fredbin_col_names, 'formats': apass.fredbin_col_types}
    data = np.asarray(container.data, dtype=dtype)

    x            = (container.rect.x_max + container.rect.x_min) / 2
    y            = (container.rect.y_max + container.rect.y_min) / 2
    zone_id      = container.zone_id
    node_id      = container.node_id
    cont_id = container.container_id

    # Iterate over the fields stored in this container and compute their average.
    field_ids  = set(data['field'])
    for field_id in field_ids:
        indexes = np.where(data['field'] == field_id)
        t_data = data[indexes]
        t_container = RectContainer(x, y, t_data, zone_id, node_id, cont_id)
        ave_data = average_container(t_container)
        output.append(ave_data)

    return output

def average_container(container):
    """Parses the measurements contained within a container and averages the data.

    If there are data in this container, this function will return a dictionary
    with the following information:

    - ra - Average Right Ascention
    - ra_sig - uncertainty on Right Ascention
    - dec - Average Declination
    - dec_sig - uncertainty on Declination
    - filter_ids - list containing the identifiers of the filters in this data set
    - mags - list containing averaged magnitudes
    - mags_sig - list containing standard deviation on magnitudes
    - num_nights - list containing the number of nights an object was observed
    - num_obs - list containing the number of observations of an object
    - container_size - pair representing size of container in (RA, DEC) order

    The lists mags, mags_sig, num_nights, and num_obs are in the same order as filter_ids.

    If there is no data in this container, the function will return an empty dictionary.
    """

    # A container should store data on precisely one star, but that data will
    # be taken through multiple photometric filters and potentially originate from
    # multiple fields (telescope pointings). In this function, we average
    # the photometry for each field and generate one output line for each field.

    output = dict()

    # Read in the data. This will be a numpy.narray object, so we can slice
    # and dice it however we would like.
    data = container.data
    dtype={'names': apass.fredbin_col_names,'formats': apass.fredbin_col_types}
    data = np.asarray(container.data, dtype=dtype)

    # add a boolean column to serve as a flag for whether or not a given datum should be used
    tmp = np.zeros(len(data))
    data = nprf.append_fields(data, ['use_data'], [tmp], dtypes=[bool])

    # if there is no data in the container or it has been moved, return an empty dictionary.
    if len(data) == 0 or container.moved == True:
        return output

    # Compute the average RA and DEC using data from all measurements
    rect = container.rect
    output['ra']             = average(data['ra'])
    output['ra_sig']         = std(data['ra'])
    output['dec']            = average(data['dec'])
    output['dec_sig']        = std(data['dec'])
    output['container_size'] = (rect.x_max - rect.x_min, rect.y_max - rect.y_min)
    output['mags']           = []
    output['mags_sig']       = []
    output['num_nights']     = []
    output['num_obs']        = []
    output['field']          = data['field'][0]

    ##
    # Filtering Stages
    ##
    data = filter_by_ccd_radius(data, sro_ccd_x_center, sro_ccd_y_center, sro_max_ccd_radius,
                                min_num_observations = sro_min_num_observations)

    # get a list of filters in numerical order
    filter_ids = sorted(set(data['filter_id']))
    output['filter_ids'] = filter_ids

    # iterate over the filters extracting the relevant statistics
    for filter_id in filter_ids:
        # extract the data in this filter
        indexes = np.where((data['filter_id'] == filter_id) & (data['use_data'] == True))
        t_data = data[indexes]

        # total number of observations, total number of nights.
        num_obs = len(t_data)
        num_nights = len(set(t_data['hjd']))

        mag     = 99.999
        mag_sig = 99.999

        if num_obs > 0:
            # magnitude and its uncertainty
            mag = average(t_data['xmag1'])
            if num_obs > 1:
                mag_sig = std(t_data['xmag1'])
            else:
                mag_sig = t_data['xerr1']

        # store the data in the output dictionary
        output['mags'].append(mag)
        output['mags_sig'].append(mag_sig)
        output['num_nights'].append(num_nights)
        output['num_obs'].append(num_obs)

    return output


def container_to_string(container_ave, photometry_format_func):
    """Converts average container dictionary to a formatted string

    Output will be as follows:
    - field_id
    - ra (deg)
    - ra_sig (arcsecond)
    - dec (deg)
    - dec_sig (in arcsecond)
    - total number of nights
    - total number of observations
    - [number of nights by filter]
    - [number of observations by filter]
    - [mag by filter]
    - [mag_sig by filter]

    returns a formatted string
    """


    # extract the relevant quantities from the dictionary
    field_id      = container_ave['field']
    ra            = container_ave['ra']
    ra_sig        = container_ave['ra_sig'] * 3600 # convert to arcsec
    dec           = container_ave['dec']
    dec_sig       = container_ave['dec_sig'] * 3600 # convert to arcsec
    total_nights  = sum(container_ave['num_nights'])
    total_num_obs = sum(container_ave['num_obs'])

    # Start with the following information
    #  # Field    RA(J2000)   raerr  DEC(J2000) decerr nobs  mobs
    #  # 0020131545  11.198035  0.477 -32.933803  0.396    3   12
    prefix = (field_id, ra, ra_sig, dec, dec_sig, total_nights, total_num_obs)
    prefix_fmt = "%25s %10.6f %6.3f %10.6f %6.3f %4i %4i "
    line = (prefix_fmt) % prefix

    # now for flags

    flag_large_mag_diff        = False
    flag_num_obs_diff          = False
    flag_large_position_errors = False
    flag_large_bounding_boxes  = False
    flags = (flag_large_mag_diff, flag_num_obs_diff, flag_large_position_errors,
             flag_large_bounding_boxes)
    flags_fmt = "%1i %1i %1i %1i "
    line += (flags_fmt) % flags

    # Next up is the photometry. We pass this off to the generation function
    line += photometry_format_func(container_ave)

    return line

def sro_photometry_format_func(container_ave):
    """Extracts the photometry from a container and returns a formatted string"""

    # now follow up with magnitudes and errors. Note that APASS splits them as
    #  (mag1, ..., magN, err1, ..., errN)

    num_obs_fmt = "%3i"
    mag_fmt = "%6.3f"

    out_mags = []
    out_mags_sig = []
    out_num_nights = []
    out_num_obs = []

    mag = 0
    mag_sig = 0
    num_nights = 0
    num_obs = 0

    # First we iterate over the number of measurements
    for filter_index in range(1, sro_num_filters + 1):

        # SRO filters are computed as follows
        # index | out <- in
        #   1 <- 3
        #   2 <- (2 - 3)
        #   3 <- 2
        #   4 <- 8
        #   5 <- 9
        #   6 <- 10
        # we mirror this here
        if filter_index == 1:
            mag, mag_sig, num_nights, num_obs = read_mags(container_ave, 3)
        elif filter_index == 2:
            mag2, mag2_sig, num_nights2, num_obs2 = read_mags(container_ave, 2)
            mag3, mag3_sig, num_nights3, num_obs3 = read_mags(container_ave, 3)

            if mag2 == 99.999 or mag3 == 99.999:
                mag = 99.999
                mag_sig = 99.999
            else:
                mag = mag2 - mag3
                mag_sig = sqrt(mag2_sig**2 + mag3_sig**2)

                num_nights = min(num_nights2, num_nights3)
                num_obs    = min(num_obs2, num_obs3)
        elif filter_index == 3:
            mag, mag_sig, num_nights, num_obs = read_mags(container_ave, 2)
        elif filter_index == 4:
            mag, mag_sig, num_nights, num_obs = read_mags(container_ave, 8)
        elif filter_index == 5:
            mag, mag_sig, num_nights, num_obs = read_mags(container_ave, 9)
        elif filter_index == 6:
            mag, mag_sig, num_nights, num_obs = read_mags(container_ave, 10)

        out_mags.append(mag)
        out_mags_sig.append(mag_sig)
        out_num_nights.append(num_nights)
        out_num_obs.append(num_obs)

    # Construct the output line
    line = ""

    # write out number of nights
    for nights in out_num_nights:
        line += " " + (num_obs_fmt) % (nights)

    # write out number of observations
    for obs in out_num_obs:
        line += " " + (num_obs_fmt) % (obs)

    # write out the magnitudes
    for mag in out_mags:
        line += " " + (mag_fmt) % (mag)

    # and the corresponding uncertainties
    for mag_sig in out_mags_sig:
        line += " " + (mag_fmt) % (mag_sig)

    return line

def read_mags(container_ave, filter_id):
    """Reads the magnitude and error from the dictionary for the specified filter."""

    filter_order = container_ave['filter_ids']
    mags = container_ave['mags']
    mags_sig = container_ave['mags_sig']
    num_nights = container_ave['num_nights']
    num_obs = container_ave['num_obs']

    mag = 99.999
    mag_sig = 99.999
    nights = 0
    obs = 0

    try:
        idx = filter_order.index(filter_id)
        mag = mags[idx]
        mag_sig = mags_sig[idx]
        nights = num_nights[idx]
        obs = num_obs[idx]
    except ValueError:
        pass

    return [mag, mag_sig, nights, obs]

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

                # keep track of which line numbers this container contributes
                line_numbers = []
                field_ids = []

                # get averages from the container. Write the results to disk while
                # also tracking the lines that this container contributes
                averages = average_by_field(container)
                for average in averages:
                    fmt_output = container_to_string(average, sro_photometry_format_func)
                    dat_file.write(fmt_output + "\n")
                    line_numbers.append(line_no)
                    field_ids.append(average['field'])
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
