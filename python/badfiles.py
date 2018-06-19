# Functions for reading exclusionary data files (a.k.a. badfiles)

# system includes
import numpy as np

# local includes
from fred import fredbin_col_names, fredbin_col_types, fredbin_col_fmt
from fred import read_fred

def read_bad_nights(filename):
    """Reads the bad nights file. Returns the result as a numpy structured
    array with key 'night_name'. The night name is formatted [n|s]XXXXX where
    'n' is for North, 's' is for South, and XXXXX is the JD night on which the
    data were acquired.

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

    return data

def read_bad_night_fields(bad_night_fields_filename):
    """Reads the bad-night-field file. Returns a numpy structued array
    with keys 'field_id', and 'night'"""

    data = read_fred(bad_night_fields_filename)

    return data[['field_id', 'night']]

