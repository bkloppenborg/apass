#!/bin/python

# Filename: dat.py
# Purpose: Functions that enable read/write operations or other interactions on .dat files.

import numpy as np

valid_formats = ['apass', 'sro']

# data format for output data.
std_dat_names = ['field_id', 'ra', 'ra_sig', 'dec', 'dec_sig',
                 'zone_id', 'node_id', 'container_id',
                 'container_width', 'container_height', 'container_area']
std_dat_types = ['S25', 'float64', 'float64', 'float64', 'float64',
                 'int32', 'int32', 'int32',
                 'float32', 'float32', 'float32']
std_dat_fmt   = ["%25s", "%10.6f", "%6.3f", "%10.6f", "%6.3f",
                 "%5i", "%5i", "%5i",
                 '%6.3f', '%6.3f', '%6.3f']

# format information specific to the APASS data
apass_num_phot     = 5
apass_obs_names    = ['num_obs_B', 'num_obs_V', 'num_obs_sg',
                    'num_obs_sr', 'num_obs_si']
apass_obs_types    = ['int32'] * apass_num_phot
apass_obs_fmt      = ['%4i'] * apass_num_phot
apass_nights_names = ['num_nights_B', 'num_nights_V', 'num_nights_sg',
                    'num_nights_sr', 'num_nights_si']
apass_nights_types = ['int32'] * apass_num_phot
apass_nights_fmt   = ['%4i'] * apass_num_phot
apass_phot_names   = ['B', 'V', 'sg', 'sr', 'si']
apass_filter_ids   = [2, 3, 8, 9, 10]
apass_phot_types   = ['float32'] * apass_num_phot
apass_phot_fmt     = ['%6.3f'] * apass_num_phot
apass_err_names    = ['B_sig', 'V_sig', 'sg_sig', 'sr_sig', 'si_sig']
apass_err_types    = ['float32'] * apass_num_phot
apass_err_fmt      = ['%6.3f'] * apass_num_phot

# format information specific to the SRO data.
sro_num_phot     = 5
sro_obs_names    = ['num_obs_B', 'num_obs_V', 'num_obs_sg',
                    'num_obs_sr', 'num_obs_si']
sro_obs_types    = ['int32'] * sro_num_phot
sro_obs_fmt      = ['%4i'] * sro_num_phot
sro_nights_names = ['num_nights_B', 'num_nights_V', 'num_nights_sg',
                    'num_nights_sr', 'num_nights_si']
sro_nights_types = ['int32'] * sro_num_phot
sro_nights_fmt   = ['%4i'] * sro_num_phot
sro_phot_names   = ['B', 'V', 'sg', 'sr', 'si']
sro_filter_ids   = [2, 3, 8, 9, 10]
sro_phot_types   = ['float32'] * sro_num_phot
sro_phot_fmt     = ['%6.3f'] * sro_num_phot
sro_err_names    = ['B_sig', 'V_sig', 'sg_sig', 'sr_sig', 'si_sig']
sro_err_types    = ['float32'] * sro_num_phot
sro_err_fmt      = ['%6.3f'] * sro_num_phot

def filter_ids(dat_type="apass"):

    filter_ids = []

    if dat_type == "apass":
        filter_ids = apass_phot_names
    elif dat_type == "sro":
        filter_ids = sro_phot_names

    return filter_ids


def filter_names(dat_type="apass"):

    filter_names = []

    if dat_type == "apass":
        filter_names = apass_phot_names
    elif dat_type == "sro":
        filter_names = sro_phot_names

    return filter_names

def select_format(dat_type = "apass"):

    # standard column formats (note, this *copies* the lists above!)
    dat_col_names = list(std_dat_names)
    dat_col_types = list(std_dat_types)
    dat_col_fmt   = list(std_dat_fmt)

    if dat_type == "apass":
        dat_col_names.extend(apass_obs_names)
        dat_col_types.extend(apass_obs_types)
        dat_col_fmt.extend(apass_obs_fmt)
        # nights numbers
        dat_col_names.extend(apass_nights_names)
        dat_col_types.extend(apass_nights_types)
        dat_col_fmt.extend(apass_nights_fmt)
        # photometry
        dat_col_names.extend(apass_phot_names)
        dat_col_types.extend(apass_phot_types)
        dat_col_fmt.extend(apass_phot_fmt)
        dat_col_names.extend(apass_err_names)
        dat_col_types.extend(apass_err_types)
        dat_col_fmt.extend(apass_err_fmt)
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
        dat_col_names.extend(sro_err_names)
        dat_col_types.extend(sro_err_types)
        dat_col_fmt.extend(sro_err_fmt)

    return dat_col_names, dat_col_types, dat_col_fmt

def read_dat(filename, dat_type="apass"):
    """Reads in a .dat file, returns it as a Numpy structured array"""

    dat_col_names, dat_col_types, dat_col_fmt = select_format(dat_type)

    """Reads in an output data file (e.g. a .dat file)"""
    dtype={'names': dat_col_names, 'formats': dat_col_types}
    data = np.loadtxt(filename, dtype=dtype)

    # convert ra/dec errors from arcmin to deg
    data['ra_sig']  /= 60
    data['dec_sig'] /= 60
    # convert container sizes to arcseconds (or square arcseconds)
    data['container_width']  /= (3600)
    data['container_height'] /= (3600)
    data['container_area']   /= (3600 * 3600)

    return data

def write_dat(filename, data, dat_type="apass"):
    """Writes a numpy structured array to disk following the specified format"""

    dat_col_names, dat_col_types, dat_col_fmt = select_format(dat_type)

    # convert RA/DEC errors to arcmin
    data['ra_sig']  *= 60
    data['dec_sig'] *= 60
    # convert container sizes to arcseconds (or square arcseconds)
    data['container_width']  *= (3600)
    data['container_height'] *= (3600)
    data['container_area']   *= (3600 * 3600)

    # write data out in the [0,360] range
    data['ra'] %= 360.0

    # compose a header
    header =  "APASS .dat file output. Format is as follows: \n"
    header += "Column names: " + ','.join(dat_col_names) + "\n"
    header += "Column types: " + ','.join(dat_col_types) + "\n"
    header += "Column formats: " + ','.join(dat_col_fmt) + "\n"

    # save to text
    np.savetxt(filename, data, fmt=dat_col_fmt, header=header)

def dicts_to_ndarray(dicts, dat_type="apass"):
    """Converts a dictionary to a structured numpy array in a .dat-friendly format"""
    output = []

    dat_col_names, dat_col_types, dat_col_fmt = select_format(dat_type)

    # convert each dict into a list in the correct order
    for d in dicts:
        # fill in any None values in the dictionary
        if(dat_type == "apass"):
            fill_none_values(d, apass_nights_names, 0)
            fill_none_values(d, apass_obs_names, 0)
            fill_none_values(d, apass_phot_names, 99.999)
            fill_none_values(d, apass_err_names, 99.999)
        if(dat_type == "sro"):
            fill_none_values(d, sro_nights_names, 0)
            fill_none_values(d, sro_obs_names, 0)
            fill_none_values(d, sro_phot_names, 99.999)
            fill_none_values(d, sro_err_names, 99.999)

        # convert the dictionary to a list in the right order
        t = []
        for k in dat_col_names:
            t.append(d[k])

        # grow the output
        output.append(tuple(t))

    # convert the output to a numpy named array
    dtype={'names': dat_col_names, 'formats': dat_col_types}
    output = np.asarray(output, dtype=dtype)
    return output

def make_dat_dict(dat_type="apass"):
    """Makes a dat-friendly dictionary of the specified dat type."""
    dat_col_names, dat_col_types, dat_col_fmt = select_format(dat_type)
    return dict.fromkeys(dat_col_names)

def fill_none_values(a_dict, keys, value):
    """Fills the keys of a_dict with the specified value"""
    for k in keys:
        if a_dict[k] == None:
            a_dict[k] = value

def list_to_ndarray(a_list, dat_type="apass"):
    """Converts a list containing data in dat column format to a numpy.ndarray"""

    # convert back to a numpy array to make the rest of the logic easier to implement
    dat_col_names, dat_col_types, dat_col_fmt = select_format(dat_type=dat_type)
    data = np.asarray(a_list, dtype={'names': dat_col_names, 'formats': dat_col_types})

    return data
