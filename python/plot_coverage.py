#!/bin/python

# system includes
import argparse
import glob
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as patches

# project includes
import dat
from quadtree import QuadTreeNode, Rect

# file-specific globals
bounds_names   = ['ra', 'dec', 'ra_min', 'ra_max', 'dec_min', 'dec_max']
bounds_types   = ['float64'] * len(bounds_names)
num_phot_names = ['num_obs_'   + s for s in dat.apass_phot_names]
num_phot_types = ['int32'] * len(num_phot_names)
num_star_names = ['num_stars_' + s for s in dat.apass_phot_names]
num_star_types = ['int32'] * len(num_star_names)
col_names = bounds_names + num_phot_names + num_star_names
col_types = bounds_types + num_phot_types + num_star_types

class SummaryLeaf(QuadTreeNode):

    def __init__(self, rect, depth, parent=None):
        QuadTreeNode.__init__(self, rect, depth, parent)

        # create a structured numpy array that will store the data.
        temp = np.zeros(len(col_names))
        dtype={'names': col_names, 'formats':col_types}
        self.data = np.asarray(temp, dtype=dtype)

    def insert(self, x, y, datum):

        for filter_name in dat.apass_phot_names:

            # construct the names of the columns we might populate
            num_obs_col    = 'num_obs_' + filter_name
            star_count_col = 'num_stars_' + filter_name

            # determine the number of observations in this filter
            num_obs = datum[num_obs_col]

            # if the number of observations is non-zero, add the results
            if num_obs > 0:
                self.data[num_obs_col] += num_obs
                self.data[star_count_col] += 1

    def get_data(self):

        data = self.data
        data['ra_min']  = self.rect.x_min
        data['ra_max']  = self.rect.x_max
        data['dec_min'] = self.rect.y_min
        data['dec_max'] = self.rect.y_max
        data['ra']      = (self.rect.x_max + self.rect.x_min) / 2
        data['dec']     = (self.rect.y_max + self.rect.y_min) / 2

        return self.data

def discrete_cmap(N, base_cmap=None):
    """Create an N-bin discrete colormap from the specified input map"""

    # Note that if base_cmap is a string or None, you can simply do
    #    return plt.cm.get_cmap(base_cmap, N)
    # The following works for string, None, or a colormap instance:

    base = plt.cm.get_cmap(base_cmap)
    color_list = base(np.linspace(0, 1, N))
    cmap_name = base.name + str(N)
    return base.from_list(cmap_name, color_list, N)

def coords_to_pixel(ra, dec, width, height):
    x = int(ra / 360 * width)
    y = int(-dec / 180 * height) + height / 2

    return [x,y]

def main():
    """Creates a coverage map which shows the number of observations in each
    node in each zone."""

    # the coverage map is (internally) represented as a (square) Numpy NDArray
    # whose dimensions are defined by the number of times the sphere is sudivided.

    max_obs = 10 # a threshold value

    parser = argparse.ArgumentParser(description='Produces a coverage map with data from all zones')
    parser.add_argument('save_dir', type=str, help="Directory in which results (dats) are saved.")
    args = parser.parse_args()

    # collect all of the dat files.
    filenames = glob.glob((args.save_dir + "/z*.dat"))
    filenames = sorted(filenames)
    if len(filenames) == 0:
        print("No DAT files found in %s" % (args.save_dir))

    # build a tree that summarizes the properties of stars in its node
    depth = 7
    bounds = Rect(0, 360, -90, 90)
    tree = QuadTreeNode(bounds, 0)
    tree.split_until(depth, leafClass=SummaryLeaf)

    # determine the size of the image we will generate
    width  = 2**depth
    height = width

    # read in the datfile entries to the tree
    for filename in filenames:
        dat_data = dat.read_dat(filename)
        print("Read %s which has %i stars" % (filename, len(dat_data)))

        for datum in dat_data:
            x = datum['ra']
            y = datum['dec']
            tree.insert(x,y, datum)

    # Extract the data from the leaves
    leaves = tree.get_leaves()
    data = [leaf.get_data() for leaf in leaves]
    data = np.concatenate(data, axis=0)

    # create a distinct number of colors for plotting
    color_map = discrete_cmap(max_obs, 'gnuplot')

    filter_names = dat.apass_phot_names
    for filter_name in filter_names:

        # set up the image
        image = np.zeros((width, height), dtype=float)

        # determine the filter and star columns
        filter_col = 'num_obs_' + filter_name
        stars_col  = 'num_stars_' + filter_name

        # extract the average number of observations for each pixel
        for datum in data:
            ra = datum['ra']
            dec = datum['dec']
            num_obs = datum[filter_col]
            num_stars = datum[stars_col]

            ave_obs = 0
            if(num_obs > 0):
                ave_obs = float(num_obs) / num_stars

            if ave_obs > 10:
                ave_obs = 10

            [x,y] = coords_to_pixel(ra, dec, width, height)

            image[y,x] = ave_obs

        # configure the plot
        fig = plt.figure(figsize=(18,9))
        im = plt.imshow(image, extent=[0, 360, -90, 90])
        plt.colorbar(im)

        #axes.set_ylim([-90, 90])
        #axes.set_xlim([0, 360])

        plt.title("APASS %s coverage" % (filter_name))
        plt.savefig("apass_" + filter_name + ".png")


if __name__ == "__main__":
    main()
