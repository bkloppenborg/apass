#!/bin/python

# system includes
import argparse
import os
import numpy as np
import numpy.lib.recfunctions as nprf
# suppress FutureWarning from np.average
import warnings
warnings.simplefilter(action = "ignore", category = FutureWarning)
# matplotlib, set to PNG output by default
#import matplotlib as mpl
#mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import NullFormatter
from scipy.stats import norm
from scipy.optimize import least_squares

# quadtree data structure
from quadtree import *
from quadtree_types import *
from apass_types import *

# apass-specific things
import apass
from make_zones import merge_polar_zones, number_zones
import dat
import comp_file

valid_ref_formats = ['sdss']

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
    for leaf in leaves:
        for container in leaf.containers:
            # skip containers with only one entry
            if len(container.data) == 1:
                continue
            else:
                comp = container.data[0]
                ref = container.data[1]

                temp = comp_file.make_comp_dict()
                temp['ra']  = ref['ra']
                temp['dec'] = ref['dec']

                try:
                    temp['field_id']     = comp['field_id']
                    temp['zone_id']      = comp['zone_id']
                    temp['node_id']      = comp['node_id']
                    temp['container_id'] = comp['container_id']
                except ValueError:
                    continue

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
    print("Wrote output to %s" % (filename))

    # 8. Generate some plots

    # first plot, residuals vs. RA
    plot_filters = ['sg', 'sr', 'si']
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
    resid =  (A * np.exp(-(x-mu)**2 / sigma)) - y
    return resid

def scatter_histogram(x, y, hist_x_max=None, hist_y_max=None,
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
    ax_scatter.scatter(x, y, marker='.')

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

def read_sdss(filename):
    """Reads in pertinent values from a SDSS-formatted file, namely """

    # d d c                     S    RAJ        DEJ          umag   e_umag gmag   e_gmag rmag   e_rmag imag   e_imag zmag   e_zmag
    # e e l SDSS9               9 Im 2000 (deg) 2000 (deg) Q (mag)  (mag)  (mag)  (mag)  (mag)  (mag)  (mag)  (mag)  (mag)  (mag) 
    # - - - ------------------- - -- ---------- ---------- - ------ ------ ------ ------ ------ ------ ------ ------ ------ ------
    # 1 + 6 J114806.59-013005.6   Im 177.027478 -01.501561 3 15.347  0.004 14.095  0.007 13.082  0.003 13.216  0.002 12.813  0.003

    # skip the first 87 and last 3 lines in the file:
    skipheader = 87
    skiptail = -3
    infile = open(filename, 'r')
    lines = infile.readlines()[skipheader:]
    lines = lines[:skiptail]

    # Iterate through the data. Extract (ra, dec, g, r, i) and any applicable uncertainties
    # Unfortunately numpy.loadtxt fails on these files due to missing fields
    # and astropy.io.ascii is confused by the header.
    data = []
    for line in lines:
        line  = line.strip()

        temp = dict()

        ra    = float(line[31:42])
        dec   = float(line[42:53])

        try:
            sg     = float(line[69:76])
            sg_sig = float(line[77:83])
        except ValueError:
            sg = 99.999
            sg_sig = 99.999

        try:
            sr     = float(line[83:90])
            sr_sig = float(line[91:97])
        except ValueError:
            sr = 99.999
            sr_sig = 99.999

        try:
            si     = float(line[97:104])
            si_sig = float(line[105:111])
        except ValueError:
            si = 99.999
            si_sig = 99.999

        data.append((ra, dec, sg, sg_sig, sr,  sr_sig, si, si_sig, 0, 0, 0))

    # return the data as a formatted list with names matching APASS fields
    col_names = ['ra', 'dec', 'sg', 'sg_sig', 'sr', 'sr_sig', 'si', 'si_sig',
                 'zone_id', 'node_id', 'container_id']
    col_types = ['float32', 'float32', 'float32', 'float32',
                 'float32', 'float32', 'float32', 'float32',
                 'int', 'int', 'int']
    dtype={'names': col_names,'formats': col_types}
    data = np.asarray(data, dtype=dtype)

    return data


if __name__ == "__main__":
    main()
