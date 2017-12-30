#!/bin/python

import numpy as np
import numpy.lib.recfunctions as nprf
from numpy import cos
from math import pi
import re
import json

# global configuration settings for the APASS Project
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
