#!/bin/python

# Filename: comp_file.py
# Purpose: Functions that enable read/write operations or other interactions
#          on .comp files And conversions to/from ndarrays.

import numpy as np

# data format for output data.
comp_std_names = ['field_id', 'ra', 'dec',
                  'zone_id', 'node_id', 'container_id']
comp_std_types = ['S25', 'float64', 'float64',
                  'int32', 'int32', 'int32']
comp_std_fmt   = ["%25s", "%10.6f", "%10.6f",
                 "%5i", "%5i", "%5i"]

comp_phot_names = ['B', 'B_V', 'V', 'sg', 'sr', 'si']
comp_phot_types = ['float32'] * 6
comp_phot_fmt   = ['%07.3f'] * 6

def select_format():

    # standard column formats (note, this *copies* the lists above!)
    names = list(comp_std_names)
    types = list(comp_std_types)
    fmt   = list(comp_std_fmt)

    # photometry
    names.extend(comp_phot_names)
    types.extend(comp_phot_types)
    fmt.extend(  comp_phot_fmt)

    return names, types, fmt

def read(filename):
    """Reads in a .dat file, returns it as a Numpy structured array"""

    col_names, col_types, col_fmt = select_format()

    """Reads in an output data file (e.g. a .dat file)"""
    dtype={'names': col_names, 'formats': col_types}
    data = np.loadtxt(filename, dtype=dtype)

    return data

def write(filename, data):
    """Writes a numpy structured array to disk following the specified format"""

    col_names, col_types, col_fmt = select_format()

    # save to text
    np.savetxt(filename, data, fmt=col_fmt)

def dicts_to_ndarray(dicts):
    """Converts a dictionary to a structured numpy array in a .dat-friendly format"""
    output = []

    col_names, col_types, col_fmt = select_format()

    # convert each dict into a list in the correct order
    for d in dicts:

        fill_none_values(d, comp_phot_names, 99.999)
        # convert the dictionary to a list in the right order
        t = []
        for k in col_names:
            t.append(d[k])

        # grow the output
        output.append(tuple(t))

    # convert the output to a numpy named array
    dtype={'names': col_names, 'formats': col_types}
    output = np.asarray(output, dtype=dtype)
    return output

def make_comp_dict():
    """Makes a dat-friendly dictionary of the specified dat type."""
    col_names, col_types, col_fmt = select_format()
    return dict.fromkeys(col_names)

def fill_none_values(a_dict, keys, value):
    """Fills the keys of a_dict with the specified value"""
    for k in keys:
        if a_dict[k] == None:
            a_dict[k] = value

def list_to_ndarray(a_list):
    """Converts a list containing data in dat column format to a numpy.ndarray"""

    # convert back to a numpy array to make the rest of the logic easier to implement
    col_names, col_types, col_fmt = select_format()
    data = np.asarray(a_list, dtype={'names': col_names, 'formats': col_types})

    return data

