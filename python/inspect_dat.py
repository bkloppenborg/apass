#!/bin/python

# System includes
import argparse
import numpy as np
# matplotlib, set to PNG output by default
import matplotlib as mpl
#mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# Local project includes
import dat

def main():

    # specify command-line arguments
    parser = argparse.ArgumentParser(
        description="Display and inspect statistical properties of .dat files.")
    parser.add_argument('format', type=str, default="apass", choices=dat.valid_formats)
    parser.add_argument('input', nargs='+')

    # parse the command line args
    args = parser.parse_args()

    data_list = []
    for filename in args.input:
        print("Reading %s" % (filename))
        data = dat.read_dat(filename, dat_type=args.format)
        data_list.extend(data.tolist())

    data = dat.list_to_ndarray(data_list, dat_type=args.format)

    # Histogram of area distributions
    #x = data['container_area'] * (3600 * 3600)
    #plt.hist(x, bins=10000)
    #plt.show()

    # Histogram of container widths
    #x = data['container_width'] * 3600
    #plt.hist(x, bins=1000)
    #plt.show()

    # Histogram of container heights
    #x = data['container_height'] * 3600
    #plt.hist(x, bins=1000)
    #plt.show()

    # Histogram of SI uncertainties for stars brighter than 16th magnitude
    #x = data['si']
    #idx = np.where(x < 16)
    #x = data['si_sig']
    #x = x[idx]
    #plt.hist(x, bins=1000)
    #plt.show()

    # Histogram of the number of observations in SI for stars brighter than 16th magnitude
    #x = data['si']
    #idx = np.where(x < 16)
    #x = data['num_obs_si']
    #bins = max(x)
    #plt.hist(x, bins=bins)
    #plt.show()


if __name__ == "__main__":
    main()
