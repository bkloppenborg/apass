#!/bin/python

import argparse
import numpy as np
import itertools
from scipy.optimize import least_squares
import matplotlib.pyplot as plt
import networkx as nx
from operator import attrgetter
import os

# for the scatter_histogram function
import matplotlib.gridspec as gridspec
from matplotlib.ticker import NullFormatter
from scipy.stats import norm

# code cleaning could remove/change this dependency
from apass import apass_save_dir

# the following are useful includes, but are not mission critical:
import time

BAD_MAG_VALUE = 90  # anything with a mag greater than 90 is bogus
MIN_MAG_SIG = 0.001

# An entry in a SRO .dat file is as follows
# prefix = field_id, ra, ra_sig, dec, dec_sig, total_nights, total_num_obs
# flags = large_mag_diff, num_obs_diff, large_position_errors, large_bounding_boxes
# num_nights = []
# num_obs = []
# mags = []
# mags_sig = []
# where [] are repeated for each of the filters.

# the SRO data format is as follows:
# field_id, ra, ra_sig, dec, dec_sig, num_nights, num_observations, V, (B-V), B, sg, sr, si
sro_filter_names =  ['V', 'B_V', 'B', 'sg', 'sr', 'si']
sro_num_filters = len(sro_filter_names)

sro_col_names = ['field_id', 'ra', 'ra_sig', 'dec', 'dec_sig', 'num_nights', 'num_observations']
sro_col_names.extend(['large_mag_diff', 'num_obs_diff',
                      'large_position_errors', 'large_bounding_boxes'])
sro_col_names.extend(['num_nights_' + s for s in sro_filter_names])
sro_col_names.extend(['num_obs_' + s for s in sro_filter_names])
sro_col_names.extend(sro_filter_names)
sro_col_names.extend([s + "_sig" for s in sro_filter_names])

sro_col_formats = ['%10s', '%10.6f', '%6.3f', '%10.6f', '%6.3f', '%4i', '%4i']
sro_col_formats.extend(['%1i'] * 4)                 # flags
sro_col_formats.extend(['%3i'] * sro_num_filters)   # num_nights
sro_col_formats.extend(['%3i'] * sro_num_filters)   # num_obs
sro_col_formats.extend(['%6.3f'] * sro_num_filters) # mags
sro_col_formats.extend(['%6.3f'] * sro_num_filters) # mags_sig

sro_col_types = ['S10', 'float64', 'float64', 'float64', 'float64', 'int', 'int']
sro_col_types.extend(['int'] * 4)                   # flags
sro_col_types.extend(['int'] * sro_num_filters)     # num_nights
sro_col_types.extend(['int'] * sro_num_filters)     # num_obs
sro_col_types.extend(['float32'] * sro_num_filters) # num_obs
sro_col_types.extend(['float32'] * sro_num_filters) # num_obs

def make_pointings(field_base_id, pointings):

    output = []
    for pointing in pointings:
        output.append('00' + str(field_base_id + pointing))

    return output

# Fields which we will merge (excluding the pointing number)
fields = []
fields.append([30800000, make_pointings(30800000, range(1,42))])
fields.append([30900000, make_pointings(30900000, range(1,42))])
fields.append([31000000, make_pointings(31000000, range(1,43))])
fields.append([31100000, make_pointings(31100000, range(1,42))])

# Matplotlib figure numbers for specific plots
fig_residuals_in = None
fig_residuals_out = None


def scatter_histogram(x, y, hist_x_max=None, hist_y_max=None,
                      xlabel="Magnitude (mag)",
                      ylabel="Residuals (stdevs)"):
    nullfmt = NullFormatter() # no labels
    bins = 200
    ylim = (-20, 20)
    xlim = (0, 20)

    plot_info = dict()

    gs = gridspec.GridSpec(3, 3)
    gs.update(wspace=0.05, hspace=0.05)
    ax_scatter = plt.subplot(gs[0:2, 0:2])
    ax_x_hist = plt.subplot(gs[2, 0:2])
    ax_y_hist = plt.subplot(gs[0:2, 2])

    ax_scatter.set_xlim(xlim)
    ax_scatter.set_ylim(ylim)
    ax_scatter.scatter(x, y, marker='.')

    data = ax_x_hist.hist(x, bins=bins, range=xlim)
    plot_info['hist_x_max'] = 1.10 * max(data[0])
    #n, bins, patches = ax_x_hist.hist(x, bins=bins, range=xlim)
    #plot_info['hist_x_max'] = 1.10 * max(n)
    #data = ax_y_hist.hist(y, bins=bins, range=ylim, orientation='horizontal')
    #plot_info['hist_y_max'] = 1.10 * max(data[0])
    n, bins, patches = ax_y_hist.hist(y, bins=bins, range=ylim, orientation='horizontal')
    plot_info['hist_y_max'] = 1.10 * max(n)

    params = [1000, 0, 1]
    x_bins = bins[0:-1]
    results = least_squares(gaussian_func, params, args=(x_bins, n))
    A, mu, sigma = results.x
    n_fit = gaussian_func(results.x, x_bins, np.zeros(len(n)))
    ax_y_hist.plot(n_fit, x_bins)
    label = '$A=%0.4f$\n$\mu=%0.4f$\n$\sigma=%0.4f$' % (A, mu, sigma)
    ax_y_hist.text(0.05, 0.95, label, transform=ax_y_hist.transAxes, fontsize=10,
        verticalalignment='top')

    # set styling properties
    # scatter plot
    ax_scatter.set_ylabel(ylabel)
    ax_scatter.grid()
    ax_scatter.xaxis.set_major_formatter(nullfmt)

    # x-histogram
    ax_x_hist.set_xlabel(xlabel)
    ax_x_hist.grid()
    ax_x_hist.set_xlim(xlim)
    if hist_x_max != None:
        ax_x_hist.set_ylim((0, hist_x_max))
    else:
        ax_x_hist.set_ylim((0, plot_info['hist_x_max']))

    # y-histogram
    ax_y_hist.set_ylim(ylim)
    if hist_y_max != None:
        ax_y_hist.set_xlim((0, hist_y_max))
    else:
        ax_y_hist.set_xlim((0, plot_info['hist_y_max']))
    ax_y_hist.yaxis.set_major_formatter(nullfmt)
    ax_y_hist.grid()
    ax_y_hist.legend()


    return plot_info

def gaussian_func(params, x, y):
    A, mu, sigma = params
    resid =  (A * np.exp(-(x-mu)**2 / sigma)) - y
    return resid


def get_positions(pointings, pointing_info):
    """Computes the positions of the vertexes in a graph given the row and
    column order of the vertexes.
    Returns a dictionary"""

    output = dict()
    for pointing in pointings:
        output[pointing] = pointing_info[pointing]

    return output

def read_sro_dat_file(filename):
    """Reads in a SRO .dat file to a numpy associative array."""
    dtype={'names': sro_col_names, 'formats': sro_col_types}
    data = np.loadtxt(filename, dtype=dtype)
    return data

def read_sro_centers(filename):
    """Reads in a SRO center text file and returns an associative array
    containing field_id -> (ra,dec) centers"""

    output = dict()
    with open(filename, 'r') as infile:
        for line in infile:
            if line[0] == "#":
                continue
            line = line.strip()
            line = line.split()

            # output[field_id] = (ra, dec)
            field_id = str(line[0].strip())
            ra       = float(line[1].strip())
            dec      = float(line[2].strip())
            output[field_id] = (ra, dec)

    return output

def main():
    global fig_residuals_in
    global fig_residuals_out

    parser = argparse.ArgumentParser(
        description="Merges SRO .dat files together into a coherent photometric reference frame.")
    parser.add_argument('input', nargs='+')
    parser.add_argument('--debug', default=False, action='store_true',
                        help="Run in debug mode")
    args = parser.parse_args()

    start = time.time()

    # read in the pointing information
    pointing_dict = read_sro_centers('data/sro_apass_centers.txt')

    # read in the data and overlaping node graph data structure
    data, G = read_data(args.input)
    # enable to see the graph corresponding to the linkages between fields
    #nx.draw_networkx(G)
    #plt.show()

    # Create a vertex that we can use to quickly look up all fields
    # belonging to a given field_base_id. Insert edges to link it all
    # together.
    for field_base_id, pointings in fields:

        field_name = get_field_name(field_base_id)
        G.add_node(field_name)
        for pointing in pointings:
            G.add_edge(field_name, pointing)

    # enable to see graph with field_base_id attachments
    #nx.draw_networkx(G)
    #plt.show()

    # now process each field and merge things together
    for field_base_id, pointings in fields:
        print "Merging field %s" % (field_base_id)

        fig_residuals_in = dict(x=list(), y=list())
        fig_residuals_out = dict(x=list(), y=list())

        field_name = get_field_name(field_base_id)
        neighbors = nx.all_neighbors(G, field_name)
        sub_G = G.subgraph(neighbors)
        # enable to see the graph
        plt.figure(num=None, figsize=(15,15))
        fig_graph_id = plt.gcf().number
        pos = get_positions(pointings, pointing_dict)
        nx.draw_networkx(sub_G, pos=pos)
        labels = nx.get_edge_attributes(sub_G,'weight')
        nx.draw_networkx_edge_labels(sub_G, pos, edge_labels=labels)
        plt.savefig(apass_save_dir + str(field_base_id) + "-graph.png")

        # find the node with the most edges and start there
        edges = sub_G.edges(data=True)
        edges = sorted(edges, key=lambda x: x[2]['weight'], reverse=True)
        start_node_id = edges[0][0]
        adj_node_id = edges[0][1]

        # merge the fields together.
        G.node[start_node_id]['merged'] = True
        print "Starting at node %s" % (start_node_id)
        # initial merging method, follows the path of greatest weight through the
        # graph.
        #merge_by_greatest_weight(start_node_id, data, sub_G)
        # second merging method, finds the neighbor that is most connected to existing
        # merged neighbors.
        #merge_by_neighbors(data, adj_node_id, G)
        # third merging method
        wavefront = list(sub_G.neighbors(start_node_id))
        merge_by_wavefront(data, wavefront, sub_G)

        # create two figures for residuals
        plt.figure()
        plt.suptitle("Output properties for field %i" % (field_base_id))
        x = fig_residuals_out['x']
        y = fig_residuals_out['y']
        plot_info = scatter_histogram(x,y)
        plt.savefig(apass_save_dir + str(field_base_id) + "-residuals_out.png")

        plt.figure()
        plt.suptitle("Input properties for field %i" % (field_base_id))
        x = fig_residuals_in['x']
        y = fig_residuals_in['y']
        scatter_histogram(x,y, hist_y_max=plot_info['hist_y_max'])
        plt.savefig(apass_save_dir + str(field_base_id) + "-residuals_in.png")


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

    global fig_residuals_in
    global fig_residuals_out

    # generate a message
    msg = ""
    for field_j in neighbors:
        msg += "%s " % (field_j)
    print " Merging %s into the frames of %s " % (field_i, msg)

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

        # don't include B-V in plots
        if filter_id == "B_V":
            continue

        x = x.tolist()

        residuals_in = delta_mag_func([0], x, y, x_sig, y_sig)
        residuals_in = residuals_in.tolist()
        fig_residuals_in['x'].extend(x)
        fig_residuals_in['y'].extend(residuals_in)

        residuals_out = results.fun
        residuals_out = residuals_out.tolist()
        fig_residuals_out['x'].extend(x)
        fig_residuals_out['y'].extend(residuals_out)


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
        edges = t_G.edges(data=True)
        for edge in edges:
            line_ids = edge[2]['line_ids']
            line_ids = [(a+num_rows, b+num_rows) for a,b in line_ids]
            edge[2]['line_ids'] = line_ids

        # merge the graphs
        merge_graphs(G, t_G)

        num_rows = len(data)

    # ensure that all nodes and edges are flagged as not merged.
    for node in G.nodes(data=True):
        node[1]['merged'] = False
    for edge in G.edges(data=True):
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

    for t_edge in t_G.edges(data=True):
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
