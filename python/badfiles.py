# Functions for reading exclusionary data files (a.k.a. badfiles)

# system includes
import numpy as np
import numpy.lib.recfunctions as nprf
import warnings

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

    # ensure that the entries are unique
    data = np.unique(data, axis=0)

    # Append a numeric 'night' column
    tmp = np.zeros(len(data))
    data = nprf.append_fields(data, ['night'], [tmp], dtypes=['int32'], usemask=False)

    # create a column containing the (numeric) night value
    for i in range(0, len(data)):
        data[i]['night'] = int(data[i]['night_name'][1:])

    # sort the data
    data = np.sort(data, order='night')

    return data

def read_bad_night_fields(bad_night_fields_filename):
    """Reads in the bad-night-field file. Returns a numpy structured
    array with the keys 'night', 'field_id' sorted in night-field order."""

    warnings.simplefilter("ignore")

    # load the data and extract a few columns
    data = read_fred(bad_night_fields_filename)
    data = data[['field_id', 'night']]

    # ensure that the entries are unique
    data = np.unique(data, axis=0)

    # sort the data
    data = np.sort(data, order=['night', 'field_id'])

    return data


