#!/bin/python

import argparse
import numpy as np
import itertools
from scipy.optimize import least_squares
import matplotlib.pyplot as plt
import networkx as nx
from operator import attrgetter

# code cleaning could remove/change this dependency
from apass import apass_save_dir

# the following are useful includes, but are not mission critical:
import time

BAD_MAG_VALUE = 99.999
MIN_MAG_SIG = 0.001

# the SRO data format is as follows:
# field_id, ra, ra_sig, dec, dec_sig, num_nights, num_observations, V, (B-V), B, sg, sr, si
sro_filter_names =  ['V', 'B_V', 'B', 'sg', 'sr', 'si']

sro_col_names = ['field_id', 'ra', 'ra_sig', 'dec', 'dec_sig', 'num_nights', 'num_observations']
sro_col_names.extend(sro_filter_names)
sro_col_names.extend([s + "_sig" for s in sro_filter_names])

sro_col_formats = ['%010i', '%10.6f', '%6.3f', '%10.6f', '%6.3f', '%4i', '%4i',
                   '%6.3f', '%6.3f', '%6.3f', '%6.3f', '%6.3f', '%6.3f',
                   '%6.3f', '%6.3f', '%6.3f', '%6.3f', '%6.3f', '%6.3f']
sro_col_types = ['int', 'float64', 'float64', 'float64', 'float64', 'int', 'int',
                  'float32', 'float32', 'float32', 'float32', 'float32', 'float32',
                  'float32', 'float32', 'float32', 'float32', 'float32', 'float32']

# Fields which we will merge (excluding the pointing number)
fields = []

# Fields 308 and 309 use the following pointing order:
row_order = [[ 1,  2,  6,  5,  9, 10],
             [ 3,  4,  8,  7, 11, 12],
             [13, 14, 18, 17, 21, 22],
             [15, 16, 20, 19, 23, 24],
             [31, 32, 36, 35, 27, 28],
             [29, 30, 34, 33, 25, 26]]
col_order = [5, 37, 7, 38, 17, 39, 19, 40, 35, 41, 33]
fields.append([30800000, row_order, col_order])
fields.append([30900000, row_order, col_order])

# whereas fields 310, 311, and 400 use this pointing order:
row_order = [[ 1,  2,  5,  6,  9, 10],
             [ 3,  4,  7,  8, 11, 12],
             [13, 14, 17, 18, 21, 22],
             [15, 16, 19, 20, 23, 24],
             [31, 32, 35, 36, 27, 28],
             [29, 30, 33, 34, 25, 26]]
col_order = [5, 37, 7, 38, 17, 39, 19, 40, 35, 41, 33]
fields.append([31000000, row_order, col_order])
fields.append([31100000, row_order, col_order])

def read_sro_dat_file(filename):
    """Reads in a SRO .dat file to a numpy associative array."""
    dtype={'names': sro_col_names, 'formats': sro_col_types}
    data = np.loadtxt(filename, dtype=dtype)
    return data

def main():

    parser = argparse.ArgumentParser(
        description="Merges SRO .dat files together into a coherent photometric reference frame.")
    parser.add_argument('input', nargs='+')
    parser.add_argument('--debug', default=False, action='store_true',
                        help="Run in debug mode")
    args = parser.parse_args()

    start = time.time()

    # read in the data and build a graph:
    data = read_data(args.input)
    G = build_graph(data, fields)

    # now process each field and merge things together
    for field_base_id, row_order, col_order in fields:
        print "Merging field %s" % (field_base_id)

        field_name = get_field_name(field_base_id)
        neighbors = nx.all_neighbors(G, field_name)
        sub_G = G.subgraph(neighbors)

        # find the node with the most edges
        edges = sub_G.edges(data=True)
        edges = sorted(edges, key=lambda x: x[2]['weight'], reverse=True)

        start_node_id = edges[0][0]
        merge_neighbors(start_node_id, data, sub_G)

        # verify that everyone was merged in
        for node_id in sub_G.nodes():
            if G.node[node_id]['merged'] != True:
                print "Warning: Node %i was not merged!" % (node_id)

    write_sro_dat('pAllZones.dat', data)

def merge_neighbors(node_id_i, data, G):

    # get a list of all edges attached to this node and sort them in descending
    # order by weight
    edges = G.edges(node_id_i, data=True)
    edges = sorted(edges, key=lambda x: x[2]['weight'], reverse=True)

    for edge in edges:
        node_id_i = edge[0]
        node_id_j = edge[1]

        # skip nodes that have already been merged
        if G.node[node_id_i]['merged'] == True \
           and G.node[node_id_j]['merged'] == True:
            continue

        # merge the pointings and mark the edge as merged
        merge_pointings(data, node_id_i, node_id_j, edge[2]['row_ids'])
        G.node[node_id_i]['merged'] = True
        G.node[node_id_j]['merged'] = True

        # recursively merge in other nodes
        merge_neighbors(node_id_j, data, G)



def merge_pointings(data, field_i, field_j, data_row_pairs):
    """ Computes a least squares fit between the data found in the rows identified
    in data_row_pairs and then merges field_j into the photometric system of
    field_i
    """

    print " Merging %i %i" % (field_i, field_j)

    # copy the data we will merge out of the main data array
    indexes = np.where(data['field_id'] == field_j)
    merge_field_data = data[indexes]

    # get the indices for the overlapping row entries
    i,j = map(list, zip(*data_row_pairs))
    for filter_id in sro_filter_names:

        x, y, x_sig, y_sig = get_good_data(data, i, j, filter_id)

        # run a least-squares fit between the two fields
        params = [0]
        results = least_squares(delta_mag_func, params, args=(x, y, x_sig, y_sig))
        params = results.x

        print "  %3s delta: %+f" % (filter_id, params[0])

        # find filter measurements in the merge field that are valid, copy them
        # out of the array, apply the delta, and copy them back.
        good_value_indexes = np.where(merge_field_data[filter_id] <= BAD_MAG_VALUE)
        good_rows = merge_field_data[good_value_indexes]
        good_rows[filter_id] += params[0]
        merge_field_data[good_value_indexes] = good_rows

    # copy the data we modified back into the main data array
    data[indexes] = merge_field_data


def get_field_name(field_id):
    return "field_" + str(field_id)

def read_data(filenames):

    # read in the data
    data = []
    for filename in filenames:
        print "Reading %s" % (filename)
        t_data = read_sro_dat_file(filename)

        # here we convert from a ndarray to list to speed up concatenation.
        data.extend(t_data.tolist())

    print "Read %i entries" % (len(data))

    # convert back to a numpy array to make the rest of the logic easier to implement
    data = np.asarray(data, dtype={'names': sro_col_names, 'formats': sro_col_types})

    return data

def build_graph(data, fields):
    """Builds a NetworkX graph data structure from the data.
    """

    # Find all overlapping entries. We can compare the RA/DEC values
    # directly because the sro-rect-to-data.py script ensures all
    # entries in the same container were written using the same RA/DEC
    # value
    ra  = data[0]['ra']
    dec = data[0]['dec']
    overlaps = [] # will store all overlaps after the loop completes
    t_overlaps = [(0, data[0]['field_id'])] # stores successive overlaps temporally
    for i in range(1, len(data)):
        if data[i]['ra'] == ra and data[i]['dec'] == dec:
            t_overlaps.append((i, data[i]['field_id']))
        else:
            for j in range(0, len(t_overlaps)):
                j_row_id, j_field_id = t_overlaps[j]
                for k in range(j + 1, len(t_overlaps)):
                    k_row_id, k_field_id = t_overlaps[k]

                    overlaps.append([j_field_id, j_row_id, k_field_id, k_row_id])

            # clear it
            t_overlaps = [(i, data[i]['field_id'])]

        # update/increment the RA?DEC
        ra  = data[i]['ra']
        dec = data[i]['dec']

    # Store the data into a ndarray to accelerate lookups later
    names = ['field_id0', 'row_id0', 'field_id1', 'row_id1']
    types = ['int', 'int', 'int', 'int']
    overlaps = np.core.records.fromrecords(overlaps, names=names, formats=types)

    # Now build the graph. Note, in the end it is possible that some nodes might
    # not be connected to any other nodes. We'll deal with that later.
    G = nx.Graph()

    # Add all of the nodes. Also create a field_name node so we can quickly
    # find unconnected nodes later
    for field_base_id, row_order, col_order in fields:

        field_name = "field_" + str(field_base_id)
        G.add_node(field_name)

        for row in row_order:
            for element in row:
                node_id = field_base_id + element
                G.add_node(node_id, merged=False)
                G.add_edge(field_name, node_id)

        for element in col_order:
            node_id = field_base_id + element
            G.add_node(node_id, merged=False)
            G.add_edge(field_name, node_id)

    # add edges to the graph
    for row in overlaps:

        src_field  = row['field_id0']
        src_row    = row['row_id0']
        dest_field = row['field_id1']
        dest_row   = row['row_id1']

        # attempt to find the edge, if it doesn't exist, create it.
        edge = None
        try:
            edge = G[src_field][dest_field]
        except:
            G.add_edge(src_field, dest_field, row_ids=[], weight=0, merged=False)
            edge = G[src_field][dest_field]

        edge['row_ids'].append((src_row, dest_row))
        edge['weight'] = len(edge['row_ids'])

    # the graph is built, we're golden.
    return G


def write_sro_dat(filename, data):
    print "Saving %s" % (filename)
    np.savetxt(filename, data, fmt=sro_col_formats)

def get_good_data(data, i, j, filter_id):

    # extract the data for the specified filter
    x     = data[i][filter_id]
    y     = data[j][filter_id]
    x_sig = data[i][filter_id + "_sig"]
    y_sig = data[j][filter_id + "_sig"]

    # remove nonsense/bad values
    indexes = np.where((x <= BAD_MAG_VALUE) & (y <= BAD_MAG_VALUE))
    x     = x[indexes]
    y     = y[indexes]
    x_sig = x_sig[indexes]
    y_sig = y_sig[indexes]

    # Enforce a minimum uncertainty
    x_sig[np.where(x_sig < MIN_MAG_SIG)] = MIN_MAG_SIG
    y_sig[np.where(y_sig < MIN_MAG_SIG)] = MIN_MAG_SIG

    return x, y, x_sig, y_sig


def delta_mag_func(params, x, y, x_sig, y_sig):
    """Objective function for computing the linear offset between two pointings"""
    delta = params[0]
    return (x - y + delta) / np.sqrt(x_sig * y_sig)

if __name__ == "__main__":
    main()
