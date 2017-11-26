#!/bin/python

import argparse
import numpy
import apass

def main():
    """Exports a .fredbin file to text"""

    parser = argparse.ArgumentParser(description='Merges .fred photometric files')
    parser.add_argument('input', nargs='+', help="Input files which will be split into zonefiles")

    args = parser.parse_args()

    for filename in args.input:
        data = apass.read_fredbin(filename)
        numpy.savetxt(filename + ".txt", data, fmt=apass.fredbin_savetxt_fmt)


if __name__ == "__main__":
    main()
