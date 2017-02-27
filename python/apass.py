#!/bin/python

import numpy as np
from numpy import cos
from math import pi
import re
import json

apass_save_dir = '/home/data/apass-test/'

# FRED files have the following format:
## STANDARD MAGNITUDES ONLY
## FILCON ver 3.0
## RA (J2000)    DEC    CCDX      CCDY     Flags      HJD   Airmass  Set     Group    Field      Filt   Mag     Error    dmag    sys night
## 102.5140180   0.2598700  1828.670    16.950 0 0 56295.536090 2.133    1         92 0020110040     8  10.7481  0.0050  0.0460    31 56295
apass_col_names = ['ra', 'dec', 'ccdx', 'ccdy', 'flag1', 'flag2', 'hjd', 'avexx', 'kset', 'group', 'star', 'filter_id', 'xmag1', 'xerr1', 'dmag', 'sys', 'night']
apass_col_types = ['float64', 'float64', 'float32', 'float32', 'bool', 'bool', 'float32', 'float32', 'int32', 'int32', 'int32', 'uint8', 'float32', 'float32', 'float32', 'int32', 'int32']

def read_fred(filename):
    dtype={'names': apass_col_names,'formats': apass_col_types}
    return np.loadtxt(filename, dtype=dtype)

def read_fredbin(filename):
    dtype={'names': apass_col_names,'formats': apass_col_types}
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

def node_from_name(filename):
    match = re.search("n[0-9]*", filename)
    if match:
        return int(match.group(0)[1:])
    return None

def container_from_name(filename):
    match = re.search("c[0-9]*", filename)
    if match:
        return int(match.group(0)[1:])
    return None

def name_zone_file(zone_id):
    return name_zone + ".fredbin"

def name_zone(zone_id):
    return "z" + str(zone_id).zfill(5)

def name_node(zone_id, node_id):
    name = "z" + str(zone_id).zfill(5) + "-" + \
           "n" + str(node_id).zfill(4)
    return name

def name_rect_file(zone_id, node_id, container_id):
    return name_rect(zone_id, node_id, container_id) + ".fredbin"

def name_rect(zone_id, node_id, container_id):
    name = "z" + str(zone_id).zfill(5) + "-" + \
           "n" + str(node_id).zfill(4) + "-" + \
           "c" + str(container_id).zfill(4)
    return name

def get_pair(list_object):
    x = float(list_object.pop(0))
    y = float(list_object.pop(0))
    return (x,y)

def make_border_info(container):
    info = dict()
    name = name_rect(container.zone_id, container.node_id, container.container_id)
    info['filename']     = name + ".fredbin"
    info['zone_id']      = container.zone_id
    info['node_id']      = container.node_id
    info['container_id'] = container.container_id
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
