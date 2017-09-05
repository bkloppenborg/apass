#!/bin/python

import argparse
import numpy as np
import itertools
from scipy.optimize import least_squares

# code cleaning could remove/change this dependency
from apass import apass_save_dir

# the following are useful includes, but are not mission critical:
import time

BAD_MAG_VALUE = 99.999

# the SRO data format is as follows:
# field_id, ra, ra_sig, dec, dec_sig, num_nights, num_observations, V, (B-V), B, sg, sr, si
sro_filter_names =  ['V', 'B_V', 'B', 'sg', 'sr', 'si']
sro_col_names = ['field_id', 'ra', 'ra_sig', 'dec', 'dec_sig', 'num_nights', 'num_observations']
sro_col_names.extend(sro_filter_names)

sro_col_formats = ['%010i', '%10.6f', '%6.3f', '%10.6f', '%6.3f', '%4i', '%4i',
           '%6.3f', '%6.3f', '%6.3f', '%6.3f', '%6.3f', '%6.3f']
sro_col_types = ['int', 'float64', 'float64', 'float64', 'float64', 'int', 'int',
                  'float32', 'float32', 'float32', 'float32', 'float32', 'float32']

# Fields which we will merge (excluding the pointing number)
fields = [30800000, 30900000, 310000000, 31100000, 40000000]
# Pointing numbers (arranged as they would appear on the sky)
pointings = [ 1,  2,  5,  6,  9, 10,
                     37,
              3,  4,  7,  9, 11, 12,
                     38,
             13, 14, 17, 18, 21, 22,
                     39,
             15, 16, 19, 20, 23, 24,
                     40,
             29, 30, 33, 34, 25, 26,
                     41,
             31, 32, 35, 36, 27, 28]
pointings = np.array(pointings)

# The order in which the pointings will be merged together.
merge_order = [(1,2), (2,5), (5,6), (6,9), (9,10),
               (5,37), (37,7),
               (7,4), (4,3), (7,9), (9,11), (11,12),
               (7,38), (38,17),
               (17,14), (14,13), (17,18), (18,21), (21,22),
               (17,39), (39,19),
               (19,16), (16,15), (19,20), (20,23), (23,24),
               (19,40), (40,33),
               (33,30), (30,29), (33,34), (34,25), (25,26),
               (33,41), (41,35),
               (35,32), (32,31), (35,36), (36,27), (27,28)]

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
    for field_base_id in fields:
        print "Processing pointings for field %i" % (field_base_id)
        # generate a list of valid field IDs:
        field_ids = field_base_id + pointings

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

        # Store the data into a ndarray
        names = ['field_id0', 'row_id0', 'field_id1', 'row_id1', 'used']
        types = ['int', 'int', 'int', 'int', 'bool']
        overlaps = np.core.records.fromrecords(overlaps, names=names, formats=types)
        print "Found %i stars that overlap between fields" % (len(overlaps))

        # merge the fields together
        for (i,j) in merge_order:
            i = i + field_base_id
            j = j + field_base_id
            merge_fields(i, j, t_data, overlaps)

        # save the results to a file
        write_sro_dat(apass_save_dir + "%i" % (field_base_id) + ".dat", t_data)
        quit()

def write_sro_dat(filename, data):
    np.savetxt(filename, data, fmt=sro_col_formats)

def merge_fields(i, j, data, overlaps):
    """Transforms field j into the photometric frame of field i using a
    least-squares fit."""

    # TODO: Pull in the uncertainties!

    print "Merging %i %i" % (i, j)

    indexes = np.where((overlaps['field_id0'] == i) & (overlaps['field_id1'] == j))
    t_overlaps = overlaps[indexes]

    for filter_id in sro_filter_names:
        data_indices_i = t_overlaps['row_id0']
        data_indices_j = t_overlaps['row_id1']

        x_vals, y_vals = get_good_data(data, data_indices_i, data_indices_j, filter_id)

        # run a least-squares fit between the two fields
        x = [0]
        results = least_squares(delta_mag_func, x, args=(x_vals, y_vals))
        x = results.x

        # apply the fit to field j
        indices = np.where(data['field_id'] == j)
        t_data = data[indices]
        t_data[filter_id] += x[0]

def get_good_data(data, i, j, filter_id):
    x = data[i][filter_id]
    y = data[j][filter_id]

    indexes = np.where((x < BAD_MAG_VALUE) & (y < BAD_MAG_VALUE))

    x = x[indexes]
    y = y[indexes]

    return x,y


def delta_mag_func(params, x, y):
    delta = params[0]
    return (x - y + delta)

if __name__ == "__main__":
    main()
