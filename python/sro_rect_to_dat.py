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

import dat

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
sro_filter_ids = [2, 3, 8, 9, 10]
sro_filter_names = ['B', 'V', 'sg', 'sr', 'si'] # however we write out V, (B-V), B, sg, sr, si)
sro_min_num_observations = 3

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
        ave_data = average_container(t_container)
        output.append(ave_data)

    return output

def average_container(container):
    """Parses the measurements contained within a container and averages the data.

    Data will be returned as a sro-formatted dictionary, created with dat.make_dat_dict

    If there is no data in this container, the function will return an empty dictionary.
    """

    # A container should store data on precisely one star, but that data will
    # be taken through multiple photometric filters and potentially originate from
    # multiple fields (telescope pointings). In this function, we average
    # the photometry for each field and generate one output line for each field.

    output = dat.make_dat_dict(dat_type="sro")

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
    output['field_id']         = data['field'][0]
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

        filter_idx = sro_filter_ids.index(filter_id)

        # Assign the filter values to the corresponding entries in the dictionary
        phot_name = sro_filter_names[filter_idx]
        phot_sig_name = phot_name + "_sig"
        obs_name = 'num_obs_' + phot_name
        night_name = 'num_nights_' + phot_name

        output[phot_name]     = float(mag)
        output[phot_sig_name] = float(mag_sig)
        output[obs_name]      = int(num_obs)
        output[night_name]    = int(num_nights)

# TODO: Compute B-V color
#    if output["B"] is not None and output["V"] is not None:
#        output["B_V"] = output["B"] - output["V"]

    return output

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

    # average the data and populate the graph with shared container information
    line_number = 0
    averages = []
    for leaf in leaves:
        for container in leaf.containers:
            # average the data
            c_aves = average_by_field(container)
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
    dat_filename = apass.apass_save_dir + "/" + zone_name + ".dat"
    averages = dat.dicts_to_ndarray(averages, dat_type = "sro")
    dat.write_dat(dat_filename, averages, dat_type="sro")

    # save the graph to a pickle file
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
