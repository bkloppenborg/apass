#!/bin/python

import argparse
import numpy
from apass import *

def main():

    parser = argparse.ArgumentParser(description='Merges .fred photometric files')
    parser.add_argument('input', nargs='+', help="Input files which will be split into zonefiles")

    args = parser.parse_args()

    for filename in args.input:
        data = read_fredbin(filename)
        numpy.savetxt(filename + ".txt", data)


if __name__ == "__main__":
    main()
