#!/bin/python

# system-level includes
import argparse
import os
from numpy import *
import glob

# APASS-specific things
import apass

# parallel processing
import multiprocessing as mp

# The output looks like this
##  Name    RA(J2000)   raerr  DEC(J2000) decerr nobs  mobs       filt  mag  err
#0020131545  11.198035  0.477 -32.933803  0.396    3   12 12.828  0.931 13.758 13.245 12.593 12.471  0.159  0.164  0.039  0.036  0.088  0.289

def process_rect(filename):

    # A rectangle file should contain data on precisely one star take in
    # multiple photmetric filters.

    # Read in the data. This will be a numpy.narray object, so we can slice
    # and dice it however we would like.
    data = apass.read_fredbin(filename)

    # RA and DEC:
    name = data['star'][0]
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

    # compute the number of observations and number of nights that made it through filtering
    num_observations = len(data)
    num_nights = len(set(data['night']))

    # construct the output string
    prefix_fmt = "%010i %10.6f %6.3f %10.6f %6.3f %4i %4i"
    mag_fmt = "%6.3f"
    # Start with the following information
    #  # Name    RA(J2000)   raerr  DEC(J2000) decerr nobs  mobs
    #  # 0020131545  11.198035  0.477 -32.933803  0.396    3   12
    output = (prefix_fmt) % (name, ra, ra_sig, dec, dec_sig, num_nights, num_observations)
    # now follow up with magnitudes and errors. Note that APASS splits them as
    #  (mag1, ..., magN, err1, ..., errN)
    for filter_id in range(0, apass.num_filters):
        mag = 99.999
        if filter_id in mags:
            mag = mags[filter_id]['mag']

        output += " " + (mag_fmt) % (mag)

    for filter_id in range(0, apass.num_filters):
        mag_sig = 99.999
        if filter_id in mags:
            mag_sig = mags[filter_id]['sig']

        output += " " + (mag_fmt) % (mag_sig)

    return output


def zone_to_data(zone_name):
    """Processes all of the rectangles found in zone. Zone should be a valid subdirectory
    of apass.apass_save_dir"""

    print "Processing zone " + zone_name

    save_dir = apass.apass_save_dir
    zone_dir = save_dir + "/" + zone_name
    rect_files = glob.glob(zone_dir + "/*.fredbin")

    # process the rectangles into data
    outfile_name = save_dir + "/" + zone_name + ".dat"
    with open(outfile_name, 'w') as outfile:
        for rectfile in rect_files:
            tmp = process_rect(rectfile)
            outfile.write(tmp + "\n")

def main():

    parser = argparse.ArgumentParser(
        description="Converts a zone's rectangle-split data into APASS photometric output data.")
    parser.add_argument('--zone', nargs='+', default=[],
                        help="Name of the zones to parse (e.g. z12345). " + \
                        "Leave blank to parse all zones.")
    parser.add_argument('-j','--jobs', type=int, help="Parallel jobs", default=1)
    parser.add_argument('--debug', type=bool, help="Run in debug mode", default=False)
    args = parser.parse_args()

    save_dir = apass.apass_save_dir

    # either parse the specified zones, or parse all zones
    zones = []
    if len(args.zone) > 0:
        zones = args.zone
    else:
        zones = [d for d in os.listdir(save_dir) if os.path.isdir(os.path.join(save_dir, d))]

    # run in debug mode
    if args.debug:
        for zone in zones:
            zone_to_data(zone)
    # run in production mode
    else:
        pool = mp.Pool(args.jobs)
        result = pool.map_async(zone_to_data, zones)
        pool.close()
        pool.join()

    #process_rect(save_dir + "z03393/z03393-n0095-c0004.fredbin")

if __name__ == "__main__":
    main()
