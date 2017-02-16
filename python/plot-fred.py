#!/bin/python

import argparse
import numpy as np

from apass import *
import matplotlib.pyplot as plt
import os.path as path

def main():

    parser = argparse.ArgumentParser(description='Plots a FRED file')
    parser.add_argument('input', nargs='+', help="Files to be plotted (individually)")

    args = parser.parse_args()

    inputs = args.input
    dtype={'names': apass_col_names,'formats': apass_col_types}

    for filename in inputs:
        print("Reading: " + filename)
        # read in the data file
        data = np.loadtxt(filename, dtype=dtype)

        print("Read " + str(data.size) + " lines")
        dec = data['dec']
        ra = data['ra'] * np.cos(dec * pi/180)

        plt.scatter(ra, dec)
        plt.show()


if __name__ == "__main__":
    main()
