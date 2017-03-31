#!/bin/python

import numpy as np
import numpy.lib.recfunctions as nprf
from numpy import cos
from math import pi
import re
import json
from copy import copy

apass_save_dir = '/home/data/apass-test/'

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

num_filters = 6

def read_fred(filename):
    dtype={'names': fred_col_names,'formats': fred_col_types}
    data = np.loadtxt(filename, dtype=dtype)

    tmp = np.zeros(len(data))
    data = nprf.append_fields(data, fredbin_extra_cols, [tmp, tmp, tmp], dtypes=fredbin_extra_types)

    return data

def read_fredbin(filename):
    dtype={'names': fredbin_col_names,'formats': fredbin_col_types}
    return np.fromfile(filename, dtype)

def get_coords_raw(ra, dec):
    ra_  = ra * cos(dec * pi / 180)
    return [ra_, dec]

def get_coords(datum):
    dec = datum['dec']
    ra_  = datum['ra'] * cos(dec * pi / 180)
    return [ra_, dec]

# store information in the following format
# write files to something like this:
#  zXXXXX.fredbin - zone file
#  zXXXXX/nYYYY-cZZZZ.fredbin - zone, node, and container

def zone_from_name(filename):
    match = re.search("z[0-9]*", filename)
    if match:
        return int(match.group(0)[1:])
    return None

def name_zone_file(zone_id):
    return name_zone(zone_id) + ".fredbin"

def name_zone_container_file(zone_id):
    return name_zone(zone_id) + '-container.fredbin'

def name_zone_border_file(zone_id):
    return name_zone(zone_id) + '-border-rects.json'

def name_zone(zone_id):
    return "z" + str(zone_id).zfill(5)

def name_zone_contrib_file(zone_id):
    return name_zone(zone_id) + "-contrib.txt"

def name_zone_json_file(zone_id):
    return name_zone(zone_id) + "-zone.json"

def name_rect(zone_id, node_id, container_id):
    output = "z" + str(zone_id).zfill(5) + \
             "n" + str(node_id).zfill(5) + \
             "c" + str(container_id).zfill(5)

    return output

def get_pair(list_object):
    x = float(list_object.pop(0))
    y = float(list_object.pop(0))
    return (x,y)

def make_border_info(container):
    info = dict()
    zone_id      = container.zone_id
    node_id      = container.node_id
    container_id = container.container_id

    name = name_rect(zone_id, node_id, container_id)

    info['name']         = name
    info['zone_id']      = zone_id
    info['node_id']      = node_id
    info['container_id'] = container_id
    info['center']       = container.rect.get_center()

    nested = {name: info}
    return nested

def save_border_info(filename, information):
    json_str = json.dumps(information)
    with open(filename, 'w') as outfile:
        outfile.write(json_str)

def load_border_info(filename):
    json_str = open(filename).read()
    info = json.loads(json_str)
    return info

def load_zone_data(tree, directory):

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

    leaves = tree.get_leaves()
    zone_id = leaves[0].zone_id
    filename = directory + '/' + name_zone_container_file(zone_id)

    # save the rectangle data to file, overwriting the contents
    with open(filename, 'a+b') as data_file:
        for leaf in leaves:
            leaf.save_data(data_file)
