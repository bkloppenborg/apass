# Functions for reading exclusionary data files (a.k.a. badfiles)

# system includes
import numpy as np

# local includes
from fred import fredbin_col_names, fredbin_col_types, fredbin_col_fmt

def read_bad_nights(filename):

    # copy the name, types, and format from the fredbin format.
    index = fredbin_col_names.index('night_name')
    names = [fredbin_col_names[index]]
    types = [fredbin_col_types[index]]
    fmt   = [fredbin_col_fmt[index]]

    dtype={'names': names,'formats': types}
    data = np.loadtxt(filename, dtype)

    return data

def read_bad_night_fields(bad_night_fields_filename):
    return []

