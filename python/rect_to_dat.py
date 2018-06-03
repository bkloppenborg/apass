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
from math import sqrt
import traceback

# APASS-specific things
import apass
from apass_types import *
from quadtree import *
from quadtree_types import *

# File I/O
import fred
import dat
import zone
import badfiles

# lockfile
import sys, os
sys.path.append(os.path.join(sys.path[0],'modules', 'FileLock', 'filelock'))
from filelock import FileLock

# parallel processing
import multiprocessing as mp
from functools import partial

# numpy array manipulation
import numpy.lib.recfunctions as nprf

# suppress FutureWarning from np.average
import warnings
warnings.simplefilter(action = "ignore", category = FutureWarning)

valid_formats = ['sro', 'apass']

# configuration settings for APASS
apass_ccd_x_center = 2048 # pixels
apass_ccd_y_center = 2048 # pixels
apass_max_ccd_radius = 2500
apass_filter_ids = dat.filter_ids(dat_type="apass")
apass_filter_names = dat.filter_names(dat_type="apass") # however we write out V, (B-V), B, sg, sr, si)
apass_min_num_observations = 3

# configuration settings for SRO
sro_ccd_x_center = 2048 # pixels
sro_ccd_y_center = 2048 # pixels
sro_max_ccd_radius = 3072.0 / 2
sro_filter_ids = dat.filter_ids(dat_type="sro")
sro_filter_names = dat.filter_names(dat_type="sro") # however we write out V, (B-V), B, sg, sr, si)
sro_min_num_observations = 3

class filter_config_data:

    def __init__(self):
        self.dat_type                  = None
        self.ccd_x_center              = None
        self.ccd_y_center              = None
        self.max_ccd_radius            = None
        self.min_num_observations      = None
        self.bad_night_filename        = None
        self.bad_night_field_filename  = None

    def load_apass_defaults(self):
        self.dat_type             = 'apass'
        self.ccd_x_center         = apass_ccd_x_center
        self.ccd_y_center         = apass_ccd_y_center
        self.max_ccd_radius       = apass_max_ccd_radius
        self.min_num_observations = apass_min_num_observations
        self.filter_ids           = apass_filter_ids
        self.filter_names         = apass_filter_names

    def load_sro_defaults(self):
        self.dat_type             = 'sro'
        self.ccd_x_center         = sro_ccd_x_center
        self.ccd_y_center         = sro_ccd_y_center
        self.max_ccd_radius       = sro_max_ccd_radius
        self.min_num_observations = sro_min_num_observations
        self.filter_ids           = sro_filter_ids
        self.filter_names         = sro_filter_names

def apply_filters(data, filter_config):

    # extract filter configuration information
    ccd_x_center              = filter_config.ccd_x_center
    ccd_y_center              = filter_config.ccd_y_center
    max_ccd_radius            = filter_config.max_ccd_radius
    min_num_observations      = filter_config.min_num_observations
    bad_night_filename        = filter_config.bad_night_filename
    bad_night_field_filename  = filter_config.bad_night_field_filename

    # load bad night data
    bad_nights        = badfiles.read_bad_nights(bad_night_filename)
    bad_nights_fields = badfiles.read_bad_night_fields(bad_night_field_filename)

    # Any given filter should ONLY set 'use_data' flags to False to avoid
    # impacting other filters.
    data = filter_bad_nights(data, bad_nights)
    data = filter_bad_night_fields(data, bad_nights_fields)
    data = filter_non_photometric_nights(data)

    data = filter_ccd_radius(data, ccd_x_center, ccd_y_center, max_ccd_radius,
                                min_num_observations = min_num_observations)

    return data

def filter_ccd_radius(data, x_center, y_center, max_ccd_radius,
                         min_num_observations = 3):
    """Sets the 'use_data' flag to False for any stars that reside OUTSIDE
    of a specified radius from the CCD center provided that a minimum number
    of observations remain. Returns the modified data array"""

    # filter out measurements outside of a specific radius from the center
    # of the CCD
    max_radius_2 = (max_ccd_radius)**2

    # Compute the (x,y) location of the stars relative to the center of the CCD
    x = data['ccdx'] - x_center
    y = data['ccdy'] - y_center
    ccd_radius_2 = x*x + y*y

    # determine measurements that reside OUTSIDE of the maximum allowable radius
    radius_indices = np.where(ccd_radius_2 > max_radius_2)

    filter_ids = set(data['filter_id'])
    for filter_id in filter_ids:
        # find measurements in this filter
        filter_indices = np.where(data['filter_id'] == filter_id)
        # find measurements outside of the maximum radius

        # determine the measurements that could be flagged bad
        bad_indices = np.intersect1d(filter_indices, radius_indices)
        num_bad = len(bad_indices)
        num_good = len(filter_indices) - num_bad

        # If the minimum number of observations remains after flagging,
        # set the 'use_data' flag to False for the bad indices.
        # Otherwise leave the data untouched.
        if num_good > min_num_observations:
            data['use_data'][bad_indices] = False

    return data

def filter_non_photometric_nights(data):
    """Sets the 'use_data' flag to False for nights that are identified as
    being non-photometric in the raw FRED fields. Returns the modified data array"""

    indexes = np.where(data['flag1'] == 1)
    data['use_data'][indexes] = False

    return data

def filter_bad_nights(data, bad_nights):
    """Sets the 'use_data' flag to False for nights that are identified as
    bad nights.
    Input:
    data - numpy array loaded using fred.read_fredbin
    bad_nights - numpy array loaded using badfiles.read_bad_nights

    Returns:
    modified numpy array in fredbin format
    """

    for i in range(0, len(bad_nights)):
        bad_night = bad_nights['night_name'][i]

        indexes = np.where(data['night_name'] == bad_night)
        data['use_data'][indexes] = False

    return data

def filter_bad_night_fields(data, bad_night_fields):
    """Sets the 'use_data' flag to False for fields on specific nights that
    have been identified as bad.  Returns the modified data array"""

    for i in range(0, len(bad_night_fields)):
        bad_night    = bad_night_fields['night_name'][i]
        bad_field_id = bad_night_fileds['field_id'][i]

        indexes = np.where((data['night_name'] == bad_night) & (data['field_id'] == bad_field_id))
        data['use_data'][indexes] = False

    return data

def average_by_field(container, filter_config):
    """Parses the measurements within a container and averages the data by field.

    This function returns a nested dictionary whose keys are the field IDs and values
    correspond to the output of average_container below.
    """

    output = []

    # Get a list of the unique field IDs in the container
    dtype={'names': fred.fredbin_col_names, 'formats': fred.fredbin_col_types}
    data = np.asarray(container.data, dtype=dtype)

    x       = (container.rect.x_max + container.rect.x_min) / 2
    y       = (container.rect.y_max + container.rect.y_min) / 2
    zone_id = container.zone_id
    node_id = container.node_id
    cont_id = container.container_id

    # Iterate over the fields stored in this container and compute their average.
    field_ids  = set(data['field'])
    for field_id in field_ids:
        indexes = np.where(data['field'] == field_id)
        t_data = data[indexes]
        t_container = RectContainer(x, y, t_data, zone_id, node_id, cont_id)
        t_container.rect = container.rect
        ave_data = average_container(t_container, filter_config)
        output.append(ave_data)

    return output

def compute_weights(data):
    """Computes the weights for APASS data. The data array should be
    a numpy array of type freddat. See fred.to_freddat for more information."""

    """Computes weights associated with APASS data"""

    num_obs = len(data)

    for i in range(0, num_obs):
        exptime   = data[i]['exposure_time']
        mag       = data[i]['xmag1']
        sig       = data[i]['xerr1']
        filter_id = data[i]['filter_id']

        weight = 1

        # here we use the raw filter numbers
        if filter_id in [7, 11, 13, 14]: # su, sz, ZS, Y
            weight = 1
        elif filter_id in [3,9,10]: # V, sr, si
            if mag < 10 and exptime >= 20:
                weight = 0
            elif mag > 10:
                if sig < 0.01:
                    weight = 3
                elif sig > 0.05:
                    weight = 1
                else:
                    weight = 2
            else:
                weight = 1
        elif filter_id in [2,8]: # B and sg
            if mag < 10 and exptime >= 40:
                weight = 0
            elif mag > 10:
                if sig < 0.01:
                    weight = 3
                elif sig > 0.05:
                    weight = 1
                else:
                    weight = 2
            else:
                weight = 1
        else:
            weight = 1

        # update the weight
        data[i]['weight'] = weight

    return data


def average_container(container, filter_config):
    """Parses the measurements contained within a container and averages the data.

    Data will be returned as a formatted dictionary, created with dat.make_dat_dict

    If there is no data in this container, the function will return an empty dictionary.
    """

    # A container should store data on precisely one star, but that data will
    # be taken through multiple photometric filters and potentially originate from
    # multiple fields (telescope pointings). In this function, we average
    # the photometry for each field and generate one output line for each field.

    output = dat.make_dat_dict(dat_type=filter_config.dat_type)

    # Read in the data. This will be a numpy.narray object, so we can slice
    # and dice it however we would like.
    data = container.data
    dtype={'names': fred.fredbin_col_names,'formats': fred.fredbin_col_types}
    data = np.asarray(container.data, dtype=dtype)

    # if there is no data in the container or it has been moved, return an empty dictionary.
    if len(data) == 0:
        return output

    # Compute the average RA and DEC using data from all measurements
    rect = container.rect
    output['field_id']         = data['field_id'][0]
    output['ra']               = average(data['ra'])
    output['ra_sig']           = std(data['ra'])
    output['dec']              = average(data['dec'])
    output['dec_sig']          = std(data['dec'])
    output['zone_id']          = data['zone_id'][0]
    output['node_id']          = data['node_id'][0]
    output['container_id']     = data['container_id'][0]
    output['container_width']  = float(rect.x_max - rect.x_min)
    output['container_height'] = float(rect.y_max - rect.y_min)
    output['container_area']   = output['container_width'] * output['container_height']

    ##
    # Filtering Stages.
    ##
    # convert the data to a freddat data format
    data = fred.to_freddat(data)

    # Apply individual filters.
    data = apply_filters(data, filter_config)
    # compute the photometric weights
    data = compute_weights(data)

    # get a list of filters in numerical order
    # NOTE: It is possible that the data contain filters not specified in
    # the standard set of filter IDs. We skip these cases below
    data_filter_ids = sorted(set(data['filter_id']))

    # iterate over the filters extracting the relevant statistics
    for filter_id in data_filter_ids:

        # skip filters not in the specified filter set
        if filter_id not in filter_config.filter_ids:
            continue

        # extract the data in this filter
        indexes = np.where((data['filter_id'] == filter_id) & (data['use_data'] == True))
        t_data = data[indexes]

        # total number of observations, total number of nights.
        num_obs = len(t_data)
        num_nights = len(set(t_data['hjd']))

        mag     = 99.999
        mag_sig = 99.999

        if num_obs >= filter_config.min_num_observations:

            # magnitude and its uncertainty
            mag = average(t_data['xmag1'], weights=t_data['weight'])
            if num_obs > 1:
                mag_sig = std(t_data['xmag1'])
            else:
                mag_sig = t_data['xerr1']

        # Assign the filter values to the corresponding entries in the dictionary
        filter_idx    = filter_config.filter_ids.index(filter_id)
        phot_name     = filter_config.filter_names[filter_idx]
        phot_sig_name = phot_name + "_sig"
        obs_name      = 'num_obs_' + phot_name
        night_name    = 'num_nights_' + phot_name

        output[phot_name]     = float(mag)
        output[phot_sig_name] = float(mag_sig)
        output[obs_name]      = int(num_obs)
        output[night_name]    = int(num_nights)

    if (output["B"] is not None and output["B"] < 99.9 and
        output["V"] is not None and output["V"] < 99.9):
        output["B_V"] = output["B"] - output["V"]
        output["B_V_sig"] = sqrt( output["B_sig"]**2 + output["V_sig"]**2 )

    return output

def sro_zone_to_dat(save_dir, filter_config, zone_container_filename):
    """Processes all of the rectangles found in zone. Zone should be a valid subdirectory
    of save_dir"""

    zone_id   = apass.zone_from_name(zone_container_filename)
    zone_name = apass.name_zone(zone_id)
    print "Processing zone " + zone_name

    # create a graph data structure for this zone
    G = nx.Graph()

    # load the zone's tree and data from disk and get the leaves
    zone_json = save_dir + apass.name_zone_json_file(zone_id)
    zone_tree = QuadTreeNode.from_file(zone_json, leafClass=RectLeaf)
    zone.load_zone_data(zone_tree, save_dir)
    leaves = zone_tree.get_leaves()

    # average the data and populate the graph with shared container information
    line_number = 0
    averages = []
    for leaf in leaves:
        for container in leaf.containers:
            # average the data
            c_aves = average_by_field(container, filter_config)
            averages.extend(c_aves)

            # populate overlapping line information
            line_numbers = []
            field_ids = []
            for c_ave in c_aves:
                line_numbers.append(line_number)
                field_ids.append(c_ave['field_id'])
                line_number += 1

            # if there are more than one field ID present, populate information
            # in the graph
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

    # write out the average information
    dat_filename = save_dir + "/" + zone_name + ".dat"
    averages = dat.dicts_to_ndarray(averages, dat_type = "sro")
    dat.write_dat(dat_filename, averages, dat_type="sro")

    # save the graph to a pickle file
    graph_filename = save_dir + "/" + zone_name + ".p"
    nx.write_gpickle(G, graph_filename)

def apass_zone_to_dat(save_dir, filter_config, zone_container_filename):
    """Process all of the rectangles found in the zone.
    With APASS data, we average on a per-container basis."""

    averages = []

    zone_id   = apass.zone_from_name(zone_container_filename)
    zone_name = apass.name_zone(zone_id)
    print "Processing zone " + zone_name

    # load the zone's tree and data from disk and get the leaves
    zone_json = save_dir + apass.name_zone_json_file(zone_id)
    zone_tree = QuadTreeNode.from_file(zone_json, leafClass=RectLeaf)
    zone.load_zone_data(zone_tree, save_dir)
    leaves = zone_tree.get_leaves()

    # average the data in each container.
    for leaf in leaves:
        for container in leaf.containers:
            # average the data
            c_ave = average_container(container, filter_config)
            averages.append(c_ave)

    # write out the average information
    dat_filename = save_dir + "/" + zone_name + ".dat"
    averages = dat.dicts_to_ndarray(averages, dat_type = "apass")
    dat.write_dat(dat_filename, averages, dat_type="apass")

def zone_to_dat(proc_func, save_dir, filter_config, zone_container_filename):
    """Wrapper function that includes exception handling and logging"""

    # input globals
    global error_filename

    try:
        proc_func(save_dir, filter_config, zone_container_filename)
    except:
        message = "ERROR: Failed to convert %s. Re-run in debug mode.\n" % (zone_container_filename)
        tb = traceback.format_exc()

        print(message)
        with FileLock(error_filename, timeout=100, delay=0.05):
            with open(error_filename, 'a') as error_file:
                error_file.write(message + "\n" + str(tb) + "\n")

def main():

    parser = argparse.ArgumentParser(
        description="Converts a containerized zone into APASS photometric output")
    parser.add_argument('-j','--jobs', type=int, help="Parallel jobs", default=4)
    parser.add_argument('format', type=str, default="apass", choices=valid_formats)
    parser.add_argument('bad_night_file', type=str,
                        help="File containing known bad nights")
    parser.add_argument('bad_night_field_file', type=str,
                        help="File containing known bad fields on specific nights")
    parser.add_argument('input', nargs='+')
    parser.add_argument('--debug', default=False, action='store_true',
                        help="Run in debug mode")

    global error_filename

    # Parse the arguments and start timing the script.
    args = parser.parse_args()
    start = time.time()

    # construct a partial function signature for execution
    save_dir = os.path.dirname(os.path.realpath(args.input[0])) + "/"

    # clear out the error file
    error_filename = save_dir + "error_rect_to_dat.txt"
    try:
        os.remove(error_filename)
    except:
        pass

    # configure the filters, assume APASS by default
    filter_config = filter_config_data()
    filter_config.load_apass_defaults()
    filter_config.bad_night_filename       = args.bad_night_file
    filter_config.bad_night_field_filename = args.bad_night_field_file

    # configure the zone-to-dat function
    ztd_func = partial(zone_to_dat, apass_zone_to_dat, save_dir, filter_config)
    if args.format == "sro":
        filter_config.load_sro_defaults()
        ztd_func = partial(zone_to_dat, sro_zone_to_dat, save_dir, filter_config)

    # run in debug mode
    if args.debug:
        for zonefile in args.input:
            ztd_func(zonefile)
    # run in production mode
    else:
        pool = mp.Pool(args.jobs)
        result = pool.imap(ztd_func, args.input)
        pool.close()
        pool.join()

    end = time.time()
    print("Time elapsed: %is" % (int(end - start)))

if __name__ == "__main__":
    main()
