#!/bin/python

import argparse
import numpy as np
import itertools
from scipy.optimize import least_squares
import matplotlib.pyplot as plt
import networkx as nx
from operator import attrgetter
import os

# code cleaning could remove/change this dependency
from apass import apass_save_dir

# the following are useful includes, but are not mission critical:
import time

BAD_MAG_VALUE = 90  # anything with a mag greater than 90 is bogus
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

def get_positions(field_base_id, row_order, col_order):
    """Computes the positions of the vertexes in a graph given the row and
    column order of the vertexes.
    Returns a dictionary"""

    num_rows = len(row_order)
    num_cols = len(row_order[0])

    pos = dict()

    for row in range(0, num_rows):
        data = row_order[row]
        for col in range(0, num_cols):
            cell_id = data[col]
            pos[field_base_id + cell_id] = (col, row)

    # find the column
    col = row_order[0].index(col_order[0])
    for row in range(0, len(col_order)):
        if row % 2 == 0:
            continue

        cell_id = col_order[row]
        pos[field_base_id + cell_id] = (col - 0.25, float(row)/2)

    return pos

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

    # read in the data and overlaping node graph data structure
    data, G = read_data(args.input)
    # enable to see the graph corresponding to the linkages between fields
    #nx.draw_networkx(G)
    #plt.show()

    # Create a vertex that we can use to quickly look up all fields
    # belonging to a given field_base_id. Insert edges to link it all
    # together.
    for field_base_id, row_order, col_order in fields:
        field_ids = []
        for row in row_order:
            for field in row:
                field_ids.append(field_base_id + field)

        for field in col_order:
            field_ids.append(field_base_id + field)

        field_ids = list(set(field_ids))
        field_name = get_field_name(field_base_id)
        G.add_node(field_name)
        for field in field_ids:
            G.add_edge(field_name, field)

    # enable to see graph with field_base_id attachments
    #nx.draw_networkx(G)
    #plt.show()

    # now process each field and merge things together
    for field_base_id, row_order, col_order in fields:
        print "Merging field %s" % (field_base_id)

        field_name = get_field_name(field_base_id)
        neighbors = nx.all_neighbors(G, field_name)
        sub_G = G.subgraph(neighbors)
        # enable to see the graph
        plt.figure(num=None, figsize=(15,15))
        pos = get_positions(field_base_id, row_order, col_order)
        nx.draw_networkx(sub_G, pos=pos)
        labels = nx.get_edge_attributes(sub_G,'weight')
        nx.draw_networkx_edge_labels(sub_G, pos, edge_labels=labels)
        plt.savefig(apass_save_dir + str(field_base_id) + ".png")

        # find the node with the most edges and start there
        edges = sub_G.edges(data=True)
        edges = sorted(edges, key=lambda x: x[2]['weight'], reverse=True)
        start_node_id = edges[0][0]
        adj_node_id = edges[0][1]

        # merge the fields together.
        G.node[start_node_id]['merged'] = True
        print "Starting at node %i" % (start_node_id)
        # initial merging method, follows the path of greatest weight through the
        # graph.
        #merge_by_greatest_weight(start_node_id, data, sub_G)
        # second merging method, finds the neighbor that is most connected to existing
        # merged neighbors.
        #merge_by_neighbors(data, adj_node_id, G)
        # third merging method
        wavefront = sub_G.neighbors(start_node_id)
        merge_by_wavefront(data, wavefront, sub_G)

        # verify that everyone was merged in
        for node_id in sub_G.nodes():
            if 'merged' in G.node[node_id] and G.node[node_id]['merged'] != True:
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
    if node_data['merged'] == True:
        return

    # find neighbors that have already been merged
    merged_neighbors = []
    other_neighbors = []
    data_row_pairs = []
    all_neighbors = nx.all_neighbors(G, field_i)
    for neighbor in all_neighbors:
        # skip non-numeric neighbor IDs (e.g. "field_XXXX")
        if isinstance(neighbor, str):
            continue

        neighbor_data = G.node[neighbor]
        if neighbor_data['merged'] == True:
            merged_neighbors.append(neighbor)
            edge_data = G.get_edge_data(field_i, neighbor)
            data_row_pairs.extend(edge_data['line_ids'])
        else:
            other_neighbors.append(neighbor)

    # merge field_i into the photometric system of its neighbors.
    merge_pointings(data, field_i, merged_neighbors, G)

    # recursively merge in all of the other neighbors
    for neighbor in other_neighbors:
        merge_by_neighbors(data, neighbor, G)

def merge_by_wavefront(data, wavefront, G):
    """
    """

    # stop merging when the wavefront is empty
    if len(wavefront) == 0:
        return

    # find the node that is most connected to merged vertices, best_node
    best_edge_count = 0
    best_node = None
    for node in wavefront:
        edge_count = 0
        for edge in G.edges(node, data=True):
            src = edge[0]
            dst = edge[1]
            weight = edge[2]['weight']

            if G.node[dst]['merged'] == True:
                edge_count += weight

        if edge_count > best_edge_count:
            best_node = node
            best_edge_count = edge_count

    # find the merged and un-merged neighbors of best_node
    all_neighbors = nx.all_neighbors(G, best_node)
    merged_neighbors = []
    unmerged_neighbors = []
    for neighbor in all_neighbors:
        if G.node[neighbor]['merged'] == True:
            merged_neighbors.append(neighbor)
        else:
            unmerged_neighbors.append(neighbor)

    # merge best_node using its merged_neighbors
    merge_pointings(data, best_node, merged_neighbors, G)

    # remove this node from the unmerged neighbors list.
    wavefront.remove(best_node)
    # append non-merged neighbors
    for neighbor in unmerged_neighbors:
        if neighbor not in wavefront:
            wavefront.append(neighbor)

    # Advance the wavefront.
    merge_by_wavefront(data, wavefront, G)

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
        if neighbor_data['merged'] == True:
            edge_data = G.get_edge_data(field_i, neighbor)
            data_row_pairs.extend(edge_data['line_ids'])

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

        print "  %3s delta: %+f, num: %i" % (filter_id, params[0], len(x))

        # find filter measurements in the merge field that are valid, copy them
        # out of the array, apply the delta, and copy them back.
        good_value_indexes = np.where(merge_field_data[filter_id] <= BAD_MAG_VALUE)
        good_rows = merge_field_data[good_value_indexes]
        good_rows[filter_id] += params[0]
        merge_field_data[good_value_indexes] = good_rows

        #plt.subplot(3,1,1)
        #point_num = np.arange(0, len(x))
        #plt.errorbar(point_num, x, xerr=x_sig, fmt='o')
        #plt.errorbar(point_num, y, xerr=y_sig, fmt='o')
        #plt.title("Input magnitude values")

        #residuals_in = delta_mag_func([0], x, y, x_sig, y_sig)
        #chi_2_in = sum(residuals_in)
        #residuals_out = results.fun
        #chi_2_out = sum(residuals_out)

        #plt.subplot(3,1,2)
        #plt.scatter(point_num, residuals_in )
        #plt.scatter(point_num, residuals_out)
        #plt.title("Residuals")

        #plt.subplot(3,1,3)
        #plt.hist(residuals_in, bins=20)
        #plt.hist(residuals_out, bins=20)
        #plt.show()

    # copy the data we modified back into the main data array
    data[indexes] = merge_field_data

    # flag the field as merged
    G.node[field_i]['merged'] = True

    return data

def get_field_name(field_id):
    return "field_" + str(field_id)

def read_data(filenames):
    """Reads in APASS/SRO formatted .dat and .p files

    Returns:
    * A numpy array containing data with the fields as described in sro_col_names
    * A networkx graph structure containing relaionships between various fields.
    """

    # read in the data
    data = []
    G = nx.Graph()
    num_rows = 0
    for filename in filenames:
        print "Reading %s" % (filename)

        base_filename = os.path.splitext(filename)[0]

        # read in the data and append it to the data list.
        t_data = read_sro_dat_file(filename)
        data.extend(t_data.tolist())

        # read in the NetworkX graph. It is formatted as follows:
        # vertex - Field ID
        # edge - dict(line_ids=[], weight=0, merged=False)
        t_G = nx.read_gpickle(base_filename + ".p")

        # update the row indicies in the edge data
        edges = t_G.edges_iter(data=True)
        for edge in edges:
            line_ids = edge[2]['line_ids']
            line_ids = [(a+num_rows, b+num_rows) for a,b in line_ids]
            edge[2]['line_ids'] = line_ids

        # merge the graphs
        merge_graphs(G, t_G)

        num_rows = len(data)

    # ensure that all nodes and edges are flagged as not merged.
    for node in G.nodes_iter(data=True):
        node[1]['merged'] = False
    for edge in G.edges_iter(data=True):
        edge[2]['merged'] = False

    print "Read %i entries" % (len(data))
    print "Nodes: %i" % (len(G.nodes()))
    print "Edges: %i" % (len(G.edges()))

    # convert back to a numpy array to make the rest of the logic easier to implement
    data = np.asarray(data, dtype={'names': sro_col_names, 'formats': sro_col_types})

    return data, G

def merge_graphs(G, t_G):
    """Merges NetworkX graph t_G into G"""
    for node in t_G.nodes():
        if node not in G.nodes():
            G.add_node(node, merged=False)

    for t_edge in t_G.edges_iter(data=True):
        t_src  = t_edge[0]
        t_dst  = t_edge[1]
        t_data = t_edge[2]

        # prohibit circular loops
        if t_src == t_dst:
            continue

        G_edge = None
        if not G.has_edge(t_src, t_dst):
            G_edge = G.add_edge(t_src, t_dst, line_ids=[], weight=0, merged=False)

        G_edge = G[t_src][t_dst]
        G_edge['weight'] += t_data['weight']
        G_edge['line_ids'].extend(t_data['line_ids'])

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
