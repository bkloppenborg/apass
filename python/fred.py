# Functions for reading and writing FRED and FREDBIN files

# system includes
import numpy as np
import numpy.lib.recfunctions as nprf
from copy import copy

# FRED files have the following format:
# STANDARD MAGNITUDES ONLY
# FILCON ver 3.3
# RA (J2000)    DEC        CCDX      CCDY  Flags   HJD      Airmass   Set      Group   Object                   Filt   Mag    Error    dmag    sys night
# 105.4134694   0.6743509  2996.030    31.010 0 0 56029.599560 1.310    1          2 10040L                        8  16.5515  0.2880  0.0391   232 56029
fred_col_names = ['ra', 'dec', 'ccdx', 'ccdy', 'flag1', 'flag2', 'hjd', 'airmass', 'set', 'group', 'field', 'filter_id', 'xmag1', 'xerr1', 'dmag', 'sys', 'night']
fred_col_types = ['float64', 'float64', 'float32', 'float32', 'bool', 'bool', 'float32', 'float32', 'int32', 'int32', 'S25', 'uint8', 'float32', 'float32', 'float32', 'int32', 'int32']

# fredbin follows the same format as fred, but also has columns for 'rect' and 'container'
fredbin_col_names = copy(fred_col_names)
fredbin_extra_col_names = ['zone_id', 'node_id', 'container_id']
fredbin_col_names.extend(fredbin_extra_col_names)
fredbin_col_types = copy(fred_col_types)
fredbin_extra_col_types = ['int32', 'int32', 'int32']
fredbin_col_types.extend(fredbin_extra_col_types)
fredbin_savetxt_fmt = ['%+011.6f', '%+011.6f', '%011.6f', '%011.6f', '%d', '%d', '%02.6f', '%02.6f', '%6i', '%6i', '%11s', '%03i', '%02.6f', '%02.6f', '%02.6f', '%6i', '%6i', '%6i', '%6i', '%6i']


def compare_fred_data(A, B):
    """Compares two FRED entries stored as numpy structured arrays.
    Returns True if they match identically, False otherwise."""

    same_point = True

    for key in fred_col_names:
        if A[key] != B[key]:
            same_point = False

    return same_point


def read_fred(filename):
    """Reads in an APASS FRED file and returns the result as a numpy structured
    array with columns as specified in fred_col_names plus fredbin.fredbin_extra_cols"""
    dtype={'names': fred_col_names,'formats': fred_col_types}
    data = np.loadtxt(filename, dtype=dtype)

    # append extra type columns for zone_id, node_id, and container_id
    tmp = np.zeros(len(data))
    data = nprf.append_fields(data, fredbin_extra_col_names, [tmp, tmp, tmp],
                              dtypes=fredbin_extra_col_types)

    return data

def read_fredbin(filename):
    """Reads in an APASS .fredbin file"""
    dtype={'names': fredbin_col_names,'formats': fredbin_col_types}
    return np.fromfile(filename, dtype)
