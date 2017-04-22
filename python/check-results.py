#!/bin/python

import argparse
import apass

import matplotlib.pyplot as plt
import numpy as np

def main():
    """Performs a zipper-style comparison between two data files."""
    parser = argparse.ArgumentParser(
        description="Converts a containerized zone into APASS photometric output")
    parser.add_argument('input', nargs='+')
    args = parser.parse_args()

    max_distance = 4.0 / (60 * 60)
    max_distance_2 = max_distance * max_distance

    fileA_name = args.input[0]
    fileB_name = args.input[1]

    fileA = apass.read_data(fileA_name)
    fileB = apass.read_data(fileB_name)

    fileA.sort(order=['ra', 'dec'])
    fileB.sort(order=['ra', 'dec'])

    print("Read in %i lines from A" % (len(fileA)))
    print("Read in %i lines from B" % (len(fileB)))

    # prepare for doing the zipper comparison
    a = 0
    b = 0
    a_max = len(fileA)
    b_max = len(fileB)
    differences = []
    while(a < a_max and b < b_max):
        # extract the the i-th (or j-th) element)
        A = fileA[a]
        B = fileB[b]
        a_ra  = A['ra']
        a_dec = A['dec']
        b_ra  = B['ra']
        b_dec = B['dec']

        distance =  (a_ra - b_ra)**2 + (a_dec - b_dec)**2
        if distance < max_distance_2:
            # the two points essentially match, advance both counters
            a += 1
            b += 1
        elif a_ra < b_ra:
            differences.append(A)
            a += 1
        elif a_ra > b_ra:
            differences.append(B)
            b += 1

    while a < a_max:
        A = file1[a]
        differences.append(A)
        a += 1

    while b < b_max:
        B = file2[b]
        differences.append(B)
        b += 1

    print("Total number of differences: %i" % (len(differences)))

    # convert the list of differences to an ndarray
    dtype={'names': apass.data_col_names, 'formats': apass.data_col_types}
    differences = np.array(differences, dtype=dtype)

    plt.scatter(fileA['ra'], fileA['dec'], color="blue", label="filstat")
    plt.scatter(fileB['ra'], fileB['dec'], color="red", label="python")
    plt.scatter(differences['ra'], differences['dec'] + 1.0 / (3600),
                color="green", marker='x', label="differences")
    plt.legend()
    plt.xlabel("RA (deg)")
    plt.ylabel("DEC (deg)")
    plt.show()

if __name__ == "__main__":
    main()
