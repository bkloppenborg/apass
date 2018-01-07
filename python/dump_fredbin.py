#!/bin/python

import argparse
import numpy

import fred

def main():
    """Exports a .fredbin file to text"""

    parser = argparse.ArgumentParser(description='Merges .fred photometric files')
    parser.add_argument('input', nargs='+', help="Input files which will be split into zonefiles")

    args = parser.parse_args()

    for filename in args.input:
        data = fred.read_fredbin(filename)
        fred.write_txt(data, filename + ".txt")

if __name__ == "__main__":
    main()
