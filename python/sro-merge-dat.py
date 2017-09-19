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
    # enable to see the graph
    #nx.draw_spring(G)
    #plt.show()

    # now process each field and merge things together
    for field_base_id, row_order, col_order in fields:
        print "Merging field %s" % (field_base_id)

        field_name = get_field_name(field_base_id)
        neighbors = nx.all_neighbors(G, field_name)
        sub_G = G.subgraph(neighbors)
        # enable to see the graph
        #nx.draw_networkx(sub_G)
        #plt.show()

        # find the node with the most edges and start there
        edges = sub_G.edges(data=True)
        edges = sorted(edges, key=lambda x: x[2]['weight'], reverse=True)
        start_node_id = edges[0][0]
        adj_node_id = edges[0][1]

        # merge the fields together.
        G.node[start_node_id]['merged'] = True
        #merge_by_greatest_weight(start_node_id, data, sub_G)
        merge_by_neighbors(data, adj_node_id, G)

        # verify that everyone was merged in
        for node_id in sub_G.nodes():
            if G.node[node_id]['merged'] != True:
                print "Warning: Node %i was not merged!" % (node_id)

        # save data specific to this field
        indexes = np.in1d(data['field_id'], sub_G.nodes())
        t_data = data[indexes]
        filename = apass_save_dir + '/p%i.dat' % (field_base_id)
        write_sro_dat(filename, t_data)

    # save data to all fields
    write_sro_dat(apass_save_dir + '/pALL.dat', data)


def merge_by_greatest_weight(node_id_i, data, G):
    """Merges neighboring nodes together by following a path of greatest wei"""

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
        merge_pointings(data, node_id_j, [node_id_i], G)

        # recursively merge in other nodes
        merge_by_greatest_weight(node_id_j, data, G)


def merge_by_neighbors(data, field_i, G):
    """Merges field_i into the photometric system of its neighboring fields
    that have already been merged."""

    # skip nodes that have already been merged
    node_data = G.node[field_i]
    if 'merged' not in node_data or node_data['merged'] == True:
        return

    # find neighbors that have already been merged
    merged_neighbors = []
    other_neighbors = []
    data_row_pairs = []
    all_neighbors = nx.all_neighbors(G, field_i)
    for neighbor in all_neighbors:
        neighbor_data = G.node[neighbor]
        if 'merged' in neighbor_data and neighbor_data['merged'] == True:
            merged_neighbors.append(neighbor)
            edge_data = G.get_edge_data(field_i, neighbor)
            data_row_pairs.extend(edge_data['row_ids'])
        else:
            other_neighbors.append(neighbor)

    # merge field_i into the photometric system of its neighbors.
    merge_pointings(data, field_i, merged_neighbors, G)

    # recursively merge in all of the other neighbors
    for neighbor in other_neighbors:
        merge_by_neighbors(data, neighbor, G)


def merge_pointings(data, field_i, neighbors, G):

    # generate a message
    msg = ""
    for field_j in neighbors:
        msg += "%i " % (field_j)
    print " Merging %i into the frames of %s " % (field_i, msg)

    # generate a list of rows from which we will extract values.
    data_row_pairs = []
    for neighbor in neighbors:
        neighbor_data = G.node[neighbor]
        if 'merged' in neighbor_data and neighbor_data['merged'] == True:
            edge_data = G.get_edge_data(field_i, neighbor)
            data_row_pairs.extend(edge_data['row_ids'])

    # select the data from field_i
    indexes = np.where(data['field_id'] == field_i)
    merge_field_data = data[indexes]

    # Extract the indices for overlapping
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

    # flag the field as merged
    G.node[field_i]['merged'] = True

    return data

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

    # Now build the graph. Note, it is possible that some nodes in which we are not
    # interested made it in, so we need to keep track of valid nodes and purge
    # invalid nodes.
    G = nx.Graph()
    valid_nodes = []

    # Add all of the nodes. Also create a field_name node so we can quickly
    # find unconnected nodes later
    for field_base_id, row_order, col_order in fields:

        field_name = "field_" + str(field_base_id)
        G.add_node(field_name)
        valid_nodes.append(field_name)

        for row in row_order:
            for element in row:
                node_id = field_base_id + element
                G.add_node(node_id, merged=False)
                G.add_edge(field_name, node_id)
                valid_nodes.append(node_id)

        for element in col_order:
            node_id = field_base_id + element
            G.add_node(node_id, merged=False)
            G.add_edge(field_name, node_id)
            valid_nodes.append(node_id)

    # make the list of valid nodes unique
    valid_nodes = list(set(valid_nodes))

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


    # remove any nodes that are not suppose to be in the graph
    for node in nx.nodes(G):
        if node not in valid_nodes:
            G.remove_node(node)

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
