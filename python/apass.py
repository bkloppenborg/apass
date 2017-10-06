#!/bin/python

import numpy as np
import numpy.lib.recfunctions as nprf
from numpy import cos
from math import pi
import re
import json
from copy import copy

# global configuration settings for the APASS Project
#apass_save_dir = '/home/data/sro-test/'
apass_save_dir = '/2/home/kloppenb/dr10-test/'
ccd_radius = 2500
min_num_observations = 3
num_filters = 6

# Number of times the sphere is subdivided
# this should be somewhat well matched to the polar_zone_cutoff below
global_depth = 6 # dRA = 5.625 dDEC = 2.8125
#global_depth = 7 # dRA = 2.8125 dDEC = 1.40625
# number of times a zone is subdivided
zone_depth = 4

polar_zone_cutoff = 88 # |dec| greater than this are considered part of the polar zone
# (static) IDs for the polar zones. Don't change these
south_zone_id = 0
north_zone_id = 1


# FRED files have the following format:
# STANDARD MAGNITUDES ONLY
# FILCON ver 3.3
# RA (J2000)    DEC        CCDX      CCDY  Flags   HJD      Airmass   Set      Group   Object                   Filt   Mag    Error    dmag    sys night
# 105.4134694   0.6743509  2996.030    31.010 0 0 56029.599560 1.310    1          2 10040L                        8  16.5515  0.2880  0.0391   232 56029
fred_col_names = ['ra', 'dec', 'ccdx', 'ccdy', 'flag1', 'flag2', 'hjd', 'airmass', 'set', 'group', 'field', 'filter_id', 'xmag1', 'xerr1', 'dmag', 'sys', 'night']
fred_col_types = ['float64', 'float64', 'float32', 'float32', 'bool', 'bool', 'float32', 'float32', 'int32', 'int32', 'S25', 'uint8', 'float32', 'float32', 'float32', 'int32', 'int32']

# fredbin follows the same format as fred, but also has columns for 'rect' and 'container'
fredbin_col_names = copy(fred_col_names)
fredbin_extra_cols = ['zone_id', 'node_id', 'container_id']
fredbin_col_names.extend(fredbin_extra_cols)
fredbin_col_types = copy(fred_col_types)
fredbin_extra_types = ['int32', 'int32', 'int32']
fredbin_col_types.extend(fredbin_extra_types)
fredbin_savetxt_fmt = ['%03.6f', '%03.6f', '%02.6f', '%02.6f', '%d', '%d', '%02.6f', '%02.6f', '%6i', '%6i', '%s25', '%i', '%02.6f', '%02.6f', '%02.6f', '%6i', '%6i', '%6i', '%6i', '%6i']

# data format for output data
data_col_names = ['name', 'ra', 'ra_err', 'dec', 'dec_err', 'nobs', 'mobs', 'mag1', 'mag2', 'mag3', 'mag4', 'mag5', 'mag6', 'mag_err1', 'mag_err2', 'mag_err3', 'mag_err4', 'mag_err5', 'mag_err6']
data_col_types = ['int', 'float64', 'float64', 'float64', 'float64', 'int', 'int', \
                  'float32', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32']

def compare_fred_data(A, B):
    """Compares two FRED entries stored as numpy structured arrays.
    Returns True if they match identically, False otherwise."""

    same_point = True

    for key in fred_col_names:
        if A[key] != B[key]:
            same_point = False

    return same_point

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
    return [ra, dec]

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

    # build a nested dictionary that gives us direct access to the
    # underlying containers. Indexing shall be
    #  container = nodes[node_id][container_id]
    node_dict = dict()
    leaves = tree.get_leaves() # these will be of type RectLeaf
    for leaf in leaves:
        node_id = leaf.node_id

        container_dict = dict()
        for container in leaf.containers:
            container_id = container.container_id

            container_dict[container_id] = container

        node_dict[node_id] = container_dict

    # load the data
    zone_id = leaves[0].zone_id
    filename = directory + '/' + name_zone_container_file(zone_id)
    data = read_fredbin(filename)

    # insert the data *directly* into the container, bipassing normal
    # restoration methods.
    for i in range(0, len(data)):
        datum = data[i]
        node_id = datum['node_id']
        container_id = datum['container_id']
        node_dict[node_id][container_id].append_data(datum)

def save_zone_data(tree, directory):
    """Saves zone data from the tree to the specified directory"""

    leaves = tree.get_leaves()
    zone_id = leaves[0].zone_id
    filename = directory + '/' + name_zone_container_file(zone_id)

    # save the rectangle data to file, overwriting the contents
    with open(filename, 'a+b') as data_file:
        for leaf in leaves:
            leaf.save_data(data_file)
