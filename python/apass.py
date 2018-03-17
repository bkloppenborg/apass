#!/bin/python

import numpy as np
import numpy.lib.recfunctions as nprf
from numpy import cos
from math import pi
import re
import json
from quadtree import Rect

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

filter_names = [
    'U','B','V','R',
    'F675W','I','HALC','F814W',
    'OPEN','Z','J','H',
    'K','KP','su','sg',
    'sr','si','sz','HA',
    'TB','TG','TR','OIII',
    'ZS','Y']

filter_ids = [
    1,2,3,4,
    4,5,4,5,
    6,3,1,2,
    3,3,7,8,
    9,10,11,12,
    24,25,26,16,
    13,14]

def filter_id_from_name(filter_name):
    """Returns the identifier number for a filter given its name."""
    idx = filter_names.index(filter_name)
    filter_id = filter_ids[idx]

    return filter_id

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


def wrap_bounds(ra, dec):
    """Wraps (ra,dec) coordinates at boundaries."""
    if dec < -90 or dec > 90:
        ra = (ra + 180) % 360
        dec = 90 - (dec + 90) % 180
    elif ra < 0 or ra > 360:
        ra = (ra + 360) % 360

    return ra, dec

def reflect_rect(rect):
    """Conditionally reflects the rectangle about the following lines:
    * RA 0 <-> 360 boundary"""

    x_min = rect.x_min
    x_max = rect.x_max
    y_min = rect.y_min
    y_max = rect.y_max

    # rectangles on the 0 <-> 360 line
    if x_min < 0:
        x_min += 360.0
        x_max += 360.0

    return Rect(x_min, x_max, y_min, y_max)

def reflect_rect_and_data(rect, data):
    """Conditionally reflects the rectangle and data bout the following lines:
    * RA 0 <-> 360 boundary"""

    # rectangles on the 0 <-> 360 line
    if rect.x_min < 0:
        for datum in data:
            datum['ra'] += 360.0

    rect = reflect_rect(rect)

    return [rect, data]
