#!/bin/python

import argparse
import numpy as np
import itertools
from scipy.optimize import least_squares
import networkx as nx

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
fields.append([40000000, row_order, col_order])

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

    # read in the data
    data = []
    for filename in args.input:
        print "Reading %s" % (filename)
        t_data = read_sro_dat_file(filename)

        # here we convert from a ndarray to list to speed up concatenation.
        data.extend(t_data.tolist())

    print "Read %i entries" % (len(data))

    # convert back to a numpy array to make the rest of the logic easier to implement
    data = np.asarray(data, dtype={'names': sro_col_names, 'formats': sro_col_types})

    # Now go through each field and optimize.
    for field_base_id, row_order, col_order in fields:
        print "Processing pointings for field %i" % (field_base_id)

        # Add the field base ID to the row order and col order lists to
        # properly reference the field + pointing ID in the data
        row_order = [[(x + field_base_id) for x in row] for row in row_order]
        col_order = [x + field_base_id for x in col_order]

        # generate a (unique) list of all of the field + pointing IDs required
        # for merging
        field_ids = list(itertools.chain(*row_order))
        field_ids.extend(col_order)
        field_ids = list(set(field_ids))

        # Select data that is only from this field (and all of its pointings)
        indexes = np.in1d(data['field_id'], field_ids)
        #print indexes
        t_data = data[indexes]

        # Find all overlapping entries. We can compare the RA/DEC values
        # directly because the sro-rect-to-data.py script ensures all
        # entries in the same container were written using the same RA/DEC
        # value
        ra  = t_data[0]['ra']
        dec = t_data[0]['dec']
        overlaps = [] # will store all overlaps after the loop completes
        t_overlaps = [(0, t_data[0]['field_id'])] # stores successive overlaps temporally
        for i in range(1, len(t_data)):
            if t_data[i]['ra'] == ra and t_data[i]['dec'] == dec:
                t_overlaps.append((i, t_data[i]['field_id']))
            else:
                for j in range(0, len(t_overlaps)):
                    j_row_id, j_field_id = t_overlaps[j]
                    for k in range(j + 1, len(t_overlaps)):
                        k_row_id, k_field_id = t_overlaps[k]

                        overlaps.append([j_field_id, j_row_id, k_field_id, k_row_id, False])

                # clear it
                t_overlaps = [(i, t_data[i]['field_id'])]

            # update/increment the RA?DEC
            ra  = t_data[i]['ra']
            dec = t_data[i]['dec']

        # Store the data into a ndarray to accelerate lookups later
        names = ['field_id0', 'row_id0', 'field_id1', 'row_id1', 'merged']
        types = ['int', 'int', 'int', 'int', 'bool']
        overlaps = np.core.records.fromrecords(overlaps, names=names, formats=types)
        print "Found %i stars that overlap between fields" % (len(overlaps))

        # enable to dump unique overlap pairs to the console
        #temp = zip(overlaps['field_id0'].tolist(), overlaps['field_id1'].tolist())
        #temp = set(temp)
        #for entry in temp:
        #    print entry

        # Merge the data together.
        # NOTE: Neither one of these functions works correctly on the 308 field.
        # This is going to be replaced, but preserved in Git.
#        merge_in_order(t_data, row_order, col_order, overlaps)
#        merge_neighboring(t_data, row_order, col_order, overlaps)

        # save the results back to the main array
        data[indexes] = t_data

        # save the results to a file
        write_sro_dat(apass_save_dir + "p%i" % (field_base_id) + ".dat", t_data)

def write_sro_dat(filename, data):
    print "Saving %s" % (filename)
    np.savetxt(filename, data, fmt=sro_col_formats)

def merge_neighboring(data, row_order, col_order, overlaps):

    start_cell = col_order[0]
    merge_neighbors(data, start_cell, overlaps)

def merge_neighbors(data, cell_id, overlaps):

    # TODO: rewrite this function to operate directly on overlaps?
    indexes = np.where(overlaps['field_id0'] == cell_id)
    merge_cells = overlaps[indexes]


    for entry in merge_cells:
        if entry['merged'] == True:
            continue

        other_id = entry['field_id1']
        merge_pointings(data, cell_id, other_id, overlaps)
        merge_neighbors(data, other_id, overlaps)

    overlaps[indexes] = merge_cells

    print overlaps
    quit()

def merge_in_order(data, row_order, col_order, overlaps):
    """Merges the field pointings together into a coherent frame in three stages:
       1) Merge together the vertical connecting strip
       2) Merge in each entry to the left of the connecting strip one at a time
       3) Merge in each entry to the right of the connecting strip one at a time.
    """
    row_width = len(row_order[0])
    num_col_entries = len(col_order)

    # merge the connecting column together first.
    for i in range(0, num_col_entries - 1):
        merge_pointings(data, col_order[i], col_order[i+1], overlaps)

    # find the index of the first column value in  the first row
    connecting_index = 0
    connecting_cell_value = col_order[0]
    for i in range(0, row_width):
        if row_order[0][i] == connecting_cell_value:
            connecting_index = i

    # merge values to the left and right in each row
    for row in row_order:

        # go to the left of the connecting column
        for i in range(connecting_index, 0, -1):
            merge_pointings(data, row[i], row[i - 1], overlaps)

        # go to the right of the connecting column
        for i in range(connecting_index, row_width - 1):
            merge_pointings(data, row[i], row[i + 1], overlaps)

    return

def merge_pointings(data, i, j, overlaps):

    print " Merging %i %i" % (i, j)

    indexes = np.where((overlaps['field_id0'] == i) & (overlaps['field_id1'] == j))
    t_overlaps = overlaps[indexes]

    for filter_id in sro_filter_names:
        data_indices_i = t_overlaps['row_id0']
        data_indices_j = t_overlaps['row_id1']

        x, y, x_sig, y_sig = get_good_data(data, data_indices_i, data_indices_j, filter_id)

        # run a least-squares fit between the two fields
        params = [0]
        results = least_squares(delta_mag_func, params, args=(x, y, x_sig, y_sig))
        params = results.x

        print "  %3s delta: %+f" % (filter_id, params[0])

        # apply the fit to field j
        indexes = np.where(data['field_id'] == j)
        t_data  = data[indexes]
        t_data[filter_id] += params[0]

        # flag it as merged
        t_overlaps['merged'] = True

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
