
#!/bin/python

# system-level includes
import argparse
import os
from numpy import *

# APASS-specific things
import apass

def process_rect(filename):

    # A rectangle file should contain data on precisely one star take in
    # multiple photmetric filters.

    # Read in the data. This will be a numpy.narray object, so we can slice
    # and dice it however we would like.
    data = apass.read_fredbin(filename)

    # RA and DEC:
    ra = average(data['ra'])
    ra_sig = std(data['ra'])
    dec = average(data['dec'])
    dec_sig = std(data['dec'])

    # filter out measurements outside of a 2048 - 100 = 1948 radius from
    # the center of the CCD. (this is for 4096x4096 images only)
    center = 2048
    radius = (2048-100)**2
    radius_2 = radius*radius

    # look up the (x,y) location relative to the center of the CCD
    x = data['ccdx'] - center
    y = data['ccdy'] - center
    ccd_radius_2 = x*x + y*y

    # filter out the measurements outside of that radius
    data = data[ccd_radius_2 < radius_2]

    # compute the magnitudes
    mags = dict()
    filters = set(data['filter_id'])
    for filter in filters:
        temp = data[data['filter_id'] == filter]
        mag = average(temp['xmag1'])
        mag_sig = std(temp['xmag1'])
        mags[filter] = {'mag': mag, 'sig': mag_sig}

    print ra, ra_sig, dec, dec_sig
    print mags


def main():

    parser = argparse.ArgumentParser(description='Merges .fred photometric files')
    args = parser.parse_args()

    save_dir = apass.apass_save_dir

    process_rect(save_dir + "z03393/z03393-n0095-c0004.fredbin")

if __name__ == "__main__":
    main()
