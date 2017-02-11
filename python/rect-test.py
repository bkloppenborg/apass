#!/bin/python

import argparse
import numpy as np

from apass import *
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import os.path as path

one_asec = 1.0 / (60*60)

class Rect():

    def __init__(self, x_min, x_max, y_min, y_max):
        self.x_min = float(x_min)
        self.x_max = float(x_max)
        self.y_min = float(y_min)
        self.y_max = float(y_max)

    def __init__(self, x_center, y_center, radius):
        self.x_min = float(x_center - radius)
        self.x_max = float(x_center + radius)
        self.y_min = float(y_center - radius)
        self.y_max = float(y_center + radius)

    # this function should be abstracted and implemented elsewhere
    def contains(self, x,y):
        if x >= self.x_min and x < self.x_max and y >= self.y_min and y < self.y_max:
            return True
        return False

    def overlaps(self, other):
        """Determines if two rectangles overlap"""
        # implemented as "when do two rectangles NOT overlap?"
        if (other.y_max < self.y_min or
            other.x_max < self.x_min or
            other.x_min > self.x_max or
            other.y_min > self.y_max):
            return False
        return True

    def merge(self, other):
        self.x_min = min(self.x_min, other.x_min)
        self.x_max = max(self.x_max, other.x_max)
        self.y_min = min(self.y_min, other.y_min)
        self.y_max = max(self.y_max, other.y_max)

    def plotrep(self):
        """Returns a patches friendly representation of self"""
        xy = (self.x_min, self.y_min)
        width = self.x_max - self.x_min
        height = self.y_max - self.y_min
        return [xy, width, height]

    def __repr__(self):
        return "[x_min: %f x_max: %f y_min: %f y_max: %f]" % (self.x_min, self.x_max, self.y_min, self.y_max)


def find_overlap(rect, rects):

    for r in rects:
        if rect.overlaps(r):
            return r

    return None

def main():

    parser = argparse.ArgumentParser(description='Plots a FRED file')
    parser.add_argument('input', nargs='+', help="Files to be plotted (individually)")

    args = parser.parse_args()

    inputs = args.input
    dtype={'names': apass_col_names,'formats': apass_col_types}

    i = 0
    for filename in inputs:
        print("Reading: " + filename)
        # read in the data file
        data = np.loadtxt(filename, dtype=dtype)

        print("Read " + str(data.size) + " lines")

        rects = []
        for datum in np.nditer(data):
            ra = datum['ra']
            dec = datum['dec']
            rect = Rect(ra, dec, one_asec)

            overlap_rect = find_overlap(rect, rects)

            if(overlap_rect is None):
                rects.append(rect)
            else:
                overlap_rect.merge(rect)

            print(i)
            i += 1


        fig = plt.scatter(data['ra'], data['dec'])

        for rect in rects:
            [xy, width, height] = rect.plotrep()
            t_rect = patches.Rectangle(xy, width, height, edgecolor='r', facecolor='none')
            fig.axes.add_patch(t_rect)

        plt.show()


if __name__ == "__main__":
    main()
