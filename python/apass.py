#!/bin/python

import numpy as np
import numpy.lib.recfunctions as nprf
from numpy import cos
from math import pi
import re
import json
from copy import copy

# global configuration settings for the APASS Project
apass_save_dir = '/home/data/apass-test/'
min_num_observations = 3

# FRED files have the following format:
## STANDARD MAGNITUDES ONLY
## FILCON ver 3.0
## RA (J2000)    DEC    CCDX      CCDY     Flags      HJD   Airmass  Set     Group    Field      Filt   Mag     Error    dmag    sys night
## 102.5140180   0.2598700  1828.670    16.950 0 0 56295.536090 2.133    1         92 0020110040     8  10.7481  0.0050  0.0460    31 56295
fred_col_names = ['ra', 'dec', 'ccdx', 'ccdy', 'flag1', 'flag2', 'hjd', 'avexx', 'kset', 'group', 'star', 'filter_id', 'xmag1', 'xerr1', 'dmag', 'sys', 'night']
fred_col_types = ['float64', 'float64', 'float32', 'float32', 'bool', 'bool', 'float32', 'float32', 'int32', 'int32', 'int32', 'uint8', 'float32', 'float32', 'float32', 'int32', 'int32']

# fredbin follows the same format as fred, but also has columns for 'rect' and 'container'
fredbin_col_names = copy(fred_col_names)
fredbin_extra_cols = ['zone_id', 'node_id', 'container_id']
fredbin_col_names.extend(fredbin_extra_cols)
fredbin_col_types = copy(fred_col_types)
fredbin_extra_types = ['int32', 'int32', 'int32']
fredbin_col_types.extend(fredbin_extra_types)

# data format for output data
data_col_names = ['name', 'ra', 'ra_err', 'dec', 'dec_err', 'nobs', 'mobs', 'mag1', 'mag2', 'mag3', 'mag4', 'mag5', 'mag6', 'mag_err1', 'mag_err2', 'mag_err3', 'mag_err4', 'mag_err5', 'mag_err6']
data_col_types = ['int', 'float64', 'float64', 'float64', 'float64', 'int', 'int', \
                  'float32', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32']

num_filters = 6

def read_data(filename):
    """Reads in an output data file (e.g. a .dat file)"""
    dtype={'names': data_col_names, 'formats': data_col_types}
    data = np.loadtxt(filename, dtype=dtype)
    return data

def read_fred(filename):
    """Reads in an APASS FRED file."""
    dtype={'names': fred_col_names,'formats': fred_col_types}
    data = np.loadtxt(filename, dtype=dtype)

    # append extra type columns for zone_id, node_id, and container_id
    tmp = np.zeros(len(data))
    data = nprf.append_fields(data, fredbin_extra_cols, [tmp, tmp, tmp], dtypes=fredbin_extra_types)

    return data

def read_fredbin(filename):
    """Reads in an APASS .fredbin file"""
    dtype={'names': fredbin_col_names,'formats': fredbin_col_types}
    return np.fromfile(filename, dtype)

def get_coords_raw(ra, dec):
    """Returns (ra', dec') = (ra*cos(dec), dec)"""
    ra_  = ra * cos(dec * pi / 180)
    return [ra_, dec]

def get_coords(datum):
    """Extracts the coordinates from a single datum and returns (ra', dec') = (ra*cos(dec), dec)"""
    dec = datum['dec']
    ra  = datum['ra']
    return get_coords_raw(ra, dec)

# store information in the following format
#  zXXXXX.fredbin - zone file

def zone_from_name(filename):
    """Extracts the zone ID from a filename"""
    match = re.search("z[0-9]*", filename)
    if match:
        return int(match.group(0)[1:])
    return None

def name_zone_file(zone_id):
    """Produce a name from a zone file given a zone ID"""
    return name_zone(zone_id) + ".fredbin"

def name_zone_container_file(zone_id):
    """Produces the name for a zone's container given a zone ID"""
    return name_zone(zone_id) + '-container.fredbin'

def name_zone_border_file(zone_id):
    """Produces the name of a zone border file given a zone ID"""
    return name_zone(zone_id) + '-border-rects.json'

def name_zone(zone_id):
    """Produces the name of a zone given a zone ID"""
    return "z" + str(zone_id).zfill(5)

def name_zone_contrib_file(zone_id):
    """Produces the name of a 'zone contribution' file given a zone ID"""
    return name_zone(zone_id) + "-contrib.txt"

def name_zone_json_file(zone_id):
    """Produces the name of a zone's JSON file given a zone ID"""
    return name_zone(zone_id) + "-zone.json"

def name_container(zone_id, node_id, container_id):
    """Produces a unique name for the container given the zone, node, and container IDs"""
    output = "z" + str(zone_id).zfill(5) + \
             "n" + str(node_id).zfill(5) + \
             "c" + str(container_id).zfill(5)

    return output

def get_pair(list_object):
    """Pops the first two elements from the list and returns them as a pair."""
    x = float(list_object.pop(0))
    y = float(list_object.pop(0))
    return (x,y)

def make_border_info(container):
    """Creates a border info entry for the specified container."""
    info = dict()
    zone_id      = container.zone_id
    node_id      = container.node_id
    container_id = container.container_id

    name = name_container(zone_id, node_id, container_id)

    info['name']         = name
    info['zone_id']      = zone_id
    info['node_id']      = node_id
    info['container_id'] = container_id
    info['center']       = container.rect.get_center()

    nested = {name: info}
    return nested

def save_border_info(filename, information):
    """Saves the border information to the specified file."""
    json_str = json.dumps(information)
    with open(filename, 'w') as outfile:
        outfile.write(json_str)

def load_border_info(filename):
    """Loads the border information from the specified file"""
    json_str = open(filename).read()
    info = json.loads(json_str)
    return info

def load_zone_data(tree, directory):
    """Restores zone data from the specified directory to the tree."""

    # build a list of all of the node IDs in the tree
    nodes = dict()
    leaves = tree.get_leaves() # these will be of type RectLeaf

    # put a reference to the leaf into a dictionary for faster indexing
    for leaf in leaves:
        node_id = leaf.node_id
        nodes[node_id] = leaf

    # read in the data and assign it to the leaves
    zone_id = leaves[0].zone_id
    filename = directory + '/' + name_zone_container_file(zone_id)
    data = read_fredbin(filename)
    for i in range(0, len(data)):
        datum = data[i]
        node_id = datum['node_id']
        nodes[node_id].insert_direct(datum)

def save_zone_data(tree, directory):
    """Saves zone data from the tree to the specified directory"""

    leaves = tree.get_leaves()
    zone_id = leaves[0].zone_id
    filename = directory + '/' + name_zone_container_file(zone_id)

    # save the rectangle data to file, overwriting the contents
    with open(filename, 'a+b') as data_file:
        for leaf in leaves:
            leaf.save_data(data_file)