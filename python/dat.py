#!/bin/python

# Filename: dat.py
# Purpose: Functions that enable read/write operations or other interactions on .dat files.

import numpy as np

# data format for output data
std_dat_names = ['field_id', 'ra', 'ra_err', 'dec', 'dec_err', 'num_obs', 'num_nights']
std_dat_types = ['S25', 'float64', 'float64', 'float64', 'float64', 'int', 'int']
std_dat_fmt   = ["%25s", "%10.6f", "%6.3f", "%10.6f", "%6.3f", "%4i", "%4i"]

# format information specific to the SRO data.
sro_num_phot     = 5
sro_obs_names    = ['n_obs_B', 'n_obs_V', 'n_obs_sg', 'n_obs_sr', 'n_obs_si']
sro_obs_types    = ['int32'] * sro_num_phot
sro_obs_fmt      = ['%4i'] * sro_num_phot
sro_nights_names = ['n_nights_B', 'n_nights_V', 'n_nights_sg', 'n_nights_sr', 'n_nights_si']
sro_nights_types = ['int32'] * sro_num_phot
sro_nights_fmt      = ['%4i'] * sro_num_phot
sro_phot_names   = ['B', 'V', 'sg', 'sr', 'si']
sro_phot_types   = ['float32'] * sro_num_phot
sro_phot_fmt     = ['%6.3f']

def select_format(dat_type = "apass"):

    # standard column formats
    dat_col_names = std_dat_names
    dat_col_types = std_dat_types
    dat_col_fmt   = std_dat_fmt

    if dat_type == "apass":
        dat_col_names.extend(apass_phot_names)
        dat_col_types.extend(apass_phot_types)
    elif dat_type == "sro":
        # observation numbers
        dat_col_names.extend(sro_obs_names)
        dat_col_types.extend(sro_obs_types)
        dat_col_fmt.extend(sro_obs_fmt)
        # nights numbers
        dat_col_names.extend(sro_nights_names)
        dat_col_types.extend(sro_nights_types)
        dat_col_fmt.extend(sro_nights_fmt)
        # photometry
        dat_col_names.extend(sro_phot_names)
        dat_col_types.extend(sro_phot_types)
        dat_col_fmt.extend(sro_phot_fmt)

    return dat_col_names, dat_col_types, dat_col_fmt

def read_dat(filename, dat_type="apass"):
    """Reads in a .dat file, returns it as a Numpy structured array"""

    dat_col_names, dat_col_types, dat_col_fmt = select_format(dat_type)

    """Reads in an output data file (e.g. a .dat file)"""
    dtype={'names': dat_col_names, 'formats': dat_col_types}
    data = np.loadtxt(filename, dtype=dtype)
    return data

def write_dat(filename, data, dat_type="apass"):
    """Writes a numpy structured array to disk following the specified format"""

    dat_col_names, dat_col_types, dat_col_fmt = select_format(dat_type)
    np.savetxt(filename, data, fmt=dat_col_fmt)

def dicts_to_ndarray(dicts, dat_type="apass"):
    """Converts a dictionary to a structured numpy array in a .dat-friendly format"""
    output = []

    dat_col_names, dat_col_types, dat_col_fmt = select_format(dat_type)

    # convert each dict into a list in the correct order
    for d in dicts:
        t = []
        for k in dat_col_names:
            t.append(d[k])
        output.append(t)

    dtype={'names': dat_col_names, 'formats': dat_col_types}
    output = np.asarray(output, dtype=dtype)
    return output

def make_dat_dict(dat_type="apass"):

    dat_col_names, dat_col_types, dat_col_fmt = select_format(dat_type)

    return dict.fromkeys(dat_col_names)



