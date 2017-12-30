#!/bin/python

# system includes
import argparse
import os
import numpy as np
# suppress FutureWarning from np.average
import warnings
warnings.simplefilter(action = "ignore", category = FutureWarning)
# matplotlib, set to PNG output by default
#import matplotlib as mpl
#mpl.use('Agg')
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import NullFormatter
from scipy.stats import norm
from scipy.optimize import least_squares
import sys

# quadtree data structure
from quadtree import *
from quadtree_types import *
from apass_types import *

# apass-specific things
import apass
from make_zones import merge_polar_zones, number_zones

# Data file I/O
import dat
import comp_file
from sdss import read_sdss
from cmc14 import read_cmc14

valid_ref_formats = ['sdss', 'cmc']

apass_filters = ['B', 'B_V', 'V', 'sg', 'sr', 'si']

def main():
    """Compares the magnitudes of two fields"""

    parser = argparse.ArgumentParser(description="Compare magnitudes in two or more files.")
    parser.add_argument('ref_format', choices = valid_ref_formats,
                        help="Reference photometric file format.")
    parser.add_argument('ref_file', help="Photometric reference file")
    parser.add_argument('comp_format', choices=dat.valid_formats,
                        help="Comparison photometric file format.")
    parser.add_argument('comp_file', nargs='+', help="Comparison photometric file(s)")

    # parse command line arguments
    args = parser.parse_args()

    # 1. create the quadtree data structure
    print("Building quadtree data structure")
    bounds = Rect(0, 360, -90, 90)
    tree = QuadTreeNode(bounds, 0)
    tree.split_until(8, leafClass=RectLeaf)
    number_zones(tree)
    merge_polar_zones(tree)
    print("Quadtree complete.")

    # 2. read in the reference file
    print("Reading reference data file")
    ref_data = None
    if args.ref_format == "sdss":
        ref_data = read_sdss(args.ref_file)
    elif args.ref_format == "cmc":
        ref_data = read_cmc14(args.ref_file)
    else:
        print("The reference input type is unsupported.")

    # 3. insert the reference data into the quadtree
    print("Inserting reference file into quadtree")
    for datum in ref_data:
        tree.insert(datum['ra'], datum['dec'], datum)

    # 4. read in the comparison file
    for filename in  args.comp_file:
        print("Processing %s" % (filename))
        comp_data = dat.read_dat(filename, dat_type=args.comp_format)

        # 5. insert the comparison file into the tree, dropping non-matching entries
        for datum in comp_data:
            tree.insert_or_drop(datum['ra'], datum['dec'], datum)

    # 6. compute the differences between the reference and comparison fields
    print("Comparing magnitudes")
    leaves = tree.get_leaves()
    results = []
    max_size = 0
    for leaf in leaves:
        for container in leaf.containers:
            # skip containers with only one entry
            if len(container.data) == 1:
                continue
            else:
                ref  = None
                comp = None

                # find the reference star. It will not have a 'field_id' column
                for star in container.data:
                    if 'field_id' not in star.dtype.names:
                        ref = star

                if ref is None:
                    continue

                # extract the reference RA/DEC
                ref_ra  = ref['ra']
                ref_dec = ref['dec']

                # find the closest star astrometrically to the reference star
                best_dist_2 = sys.float_info.max
                for star in container.data:
                    if star == ref:
                        continue

                    if 'field_id' not in star.dtype.names:
                        continue

                    star_ra  = star['ra']
                    star_dec = star['dec']
                    dist_2   = (ref_ra - star_ra)**2 + (ref_dec - star_dec)**2
                    if dist_2 < best_dist_2:
                        best_dist_2 = dist_2
                        comp = star

                if comp is None:
                    continue

                # build up the comparison dictionary.
                temp = comp_file.make_comp_dict()
                temp['ra']           = ref_ra
                temp['dec']          = ref_dec
                temp['field_id']     = comp['field_id']
                temp['zone_id']      = comp['zone_id']
                temp['node_id']      = comp['node_id']
                temp['container_id'] = comp['container_id']

                # compute the deltas for each magnitude
                # A data set might not have some data, so try-catch it
                for key in apass_filters:
                    delta = 99.999
                    try:
                        delta = ref[key] - comp[key]
                    except ValueError:
                        pass

                    temp[key] = delta

                # append this calculation to the results
                results.append(temp)

    # 7. write out the results:
    basename, extension = os.path.splitext(args.ref_file)
    filename = basename + ".comp"
    results = comp_file.dicts_to_ndarray(results)
    comp_file.write(filename, results)
    print("Wrote %i records to to %s" % (len(results), filename))

    # 8. Generate some plots
    plot_filters = []
    if args.ref_format == "sdss":
        plot_filters = ['sg', 'sr', 'si']
    elif args.ref_format == "cmc":
        plot_filters = ['sr']

    field_ids = set(results['field_id'])
    colors = np.zeros(len(results))
    color = 0
    for field_id in field_ids:
        indexes = np.where((results['field_id'] == field_id))
        colors[indexes] = color
        color += 1

    #for filter_name in plot_filters:
    #    x = results['ra']
    #    y = results['dec']
    #    z = results[filter_name]
    #    fig = plt.figure()
    #    ax = fig.add_subplot(111, projection='3d')
    #    ax.scatter(x,y,z, c=colors)
    #    ax.set_zlim([-1,1])
    #    plt.show()

    # first plot, residuals vs. RA
    x = results['ra']
    for filter_name in plot_filters:
        y = results[filter_name]
        scatter_histogram(x, y, xlabel="RA", ylabel=filter_name + " residuals", ylim=(-1,1))
        title = filter_name + '-vs-RA'
        filename = basename + '-' + title + '.png'
        plt.suptitle(title)
        plt.savefig(filename)

    # now plot residuals vs. DEC
    x = results['dec']
    for filter_name in plot_filters:
        y = results[filter_name]
        scatter_histogram(x, y, xlabel="DEC", ylabel=filter_name + " residuals", ylim=(-1,1))
        title = filter_name + '-vs-DEC'
        filename = basename + '-' + title + '.png'
        plt.suptitle(title)
        plt.savefig(filename)


def gaussian_func(params, x, y):
    A, mu, sigma = abs(params)
    resid = (A * np.exp(-(x-mu)**2 / sigma)) - y
    return resid

def linear_func(params, x, y):
    m, b = params
    resid = y - (m*x + b)
    return resid

def scatter_histogram(x, y, hist_x_max=None, hist_y_max=None,
                      colors="#1f77b4",
                      xlabel="Magnitude (mag)",
                      ylabel="Residuals (stdevs)",
                      ylim = (-20, 20)):
    nullfmt = NullFormatter() # no labels
    bins = 200

    plot_info = dict()

    gs = gridspec.GridSpec(1, 3)
    gs.update(wspace=0.05, hspace=0.05)
    ax_scatter = plt.subplot(gs[0:2, 0:2])
    ax_y_hist = plt.subplot(gs[0:2, 2])

    #ax_scatter.set_xlim(xlim)
    ax_scatter.set_ylim(ylim)
    ax_scatter.scatter(x, y, c=colors, marker='.')

    # plot a line using only "good" values
    indexes = np.where((-0.5 < y) & (y < 0.5))
    line_x = x[indexes]
    line_y = y[indexes]
    params = [0, 0]
    result = least_squares(linear_func, params,  args=(line_x, line_y))
    m,b = result.x
    line_y = m * line_x + b
    ax_scatter.plot(line_x, line_y, color='black')
    delta_mag = (max(line_x) - min(line_x)) * m
    label = '$m=%0.4f$\n$\Delta mag = %0.4f$' % (m, delta_mag)
    ax_scatter.text(0.05, 0.95, label,
                    transform=ax_scatter.transAxes,
                    fontsize=10, verticalalignment='top')

    n, bins, patches = ax_y_hist.hist(y, bins=bins, range=ylim, orientation='horizontal')
    plot_info['hist_y_max'] = 1.10 * max(n)

    params = [1000, 0, 1]
    x_bins = bins[0:-1]
    results = least_squares(gaussian_func, params, args=(x_bins, n))
    A, mu, sigma = results.x
    n_fit = gaussian_func(results.x, x_bins, np.zeros(len(n)))
    ax_y_hist.plot(n_fit, x_bins)
    label = '$A=%0.4f$\n$\mu=%0.4f$\n$\sigma=%0.4f$' % (A, mu, sigma)
    ax_y_hist.text(0.05, 0.95, label, transform=ax_y_hist.transAxes, fontsize=10,
        verticalalignment='top')

    # set styling properties
    # scatter plot
    ax_scatter.set_ylabel(ylabel)
    ax_scatter.grid()
    ax_scatter.set_xlabel(xlabel)

    # y-histogram
    ax_y_hist.set_ylim(ylim)
    if hist_y_max != None:
        ax_y_hist.set_xlim((0, hist_y_max))
    else:
        ax_y_hist.set_xlim((0, plot_info['hist_y_max']))
    ax_y_hist.yaxis.set_major_formatter(nullfmt)
    ax_y_hist.grid()
    ax_y_hist.legend()


    return plot_info

if __name__ == "__main__":
    main()
