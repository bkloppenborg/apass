# Functions for reading exclusionary data files (a.k.a. badfiles)

# system includes
import numpy as np

# local includes
from fred import fredbin_col_names, fredbin_col_types, fredbin_col_fmt
from fred import read_fred

def read_bad_nights(filename):
    """Reads in the bad nights file. Returns the result as a SORTED numpy
    structured array with the key 'night_name'

    The 'night_name' column is formatted as [n|s]XXXXX where 'n' is for North,
    's' is for South, and XXXXX is the JD night on which the data were
    acquired.

    Within the APASS pipeline, the night name should correspond to the FRED
    filename, but this is dependent upon naming conventions outside of this
    pipeline."""

    # copy the name, types, and format from the fredbin format.
    index = fredbin_col_names.index('night_name')
    names = [fredbin_col_names[index]]
    types = [fredbin_col_types[index]]
    fmt   = [fredbin_col_fmt[index]]

    dtype={'names': names,'formats': types}
    data = np.loadtxt(filename, dtype)

    # sort the data
    data = np.sort(data, order='night_name')

    return data

def read_bad_night_fields(bad_night_fields_filename):
    """Reads in the bad-night-field file. Returns a numpy structured
    array with the keys 'night', 'field_id' sorted in night-field order."""

    # load the data and extract a few columns
    data = read_fred(bad_night_fields_filename)
    data = data[['field_id', 'night']]

    # sort the data
    data = np.sort(data, order=['night', 'field_id'])

    return data


