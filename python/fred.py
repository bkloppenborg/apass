# Functions for reading and writing FRED and FREDBIN files

# system includes
import os
import numpy as np
import numpy.lib.recfunctions as nprf
from copy import copy

# FRED files have the following format:
# STANDARD MAGNITUDES ONLY
# FILCON ver 3.3
# RA (J2000)    DEC        CCDX      CCDY  Flags   HJD      Airmass   Set      Group   Object                   Filt   Mag    Error    dmag    sys night
# 105.4134694   0.6743509  2996.030    31.010 0 0 56029.599560 1.310    1          2 10040L                        8  16.5515  0.2880  0.0391   232 56029
fred_col_names = ['ra', 'dec', 'ccdx', 'ccdy', 'flag1', 'flag2', 'hjd', 'airmass', 'set', 'group', 'field_id', 'filter_id', 'xmag1', 'xerr1', 'dmag', 'sys', 'night', 'exposure_time']
fred_col_types = ['float64', 'float64', 'float32', 'float32', 'bool', 'bool', 'float32', 'float32', 'int32', 'int32', 'S25', 'uint8', 'float32', 'float32', 'float32', 'int32', 'int32', 'float32']

# fredbin follows the same format as fred, but also has columns for 'rect' and 'container'
fredbin_col_names = copy(fred_col_names)
fredbin_extra_col_names = ['zone_id', 'node_id', 'container_id', 'night_name']
fredbin_col_names.extend(fredbin_extra_col_names)
fredbin_col_types = copy(fred_col_types)
fredbin_extra_col_types = ['int32', 'int32', 'int32', 'S7']
fredbin_col_types.extend(fredbin_extra_col_types)
fredbin_col_fmt = ['%+011.6f', '%+011.6f', '%011.6f', '%011.6f', '%d', '%d', '%02.6f', '%02.6f', '%6i', '%6i', '%11s', '%03i', '%02.6f', '%02.6f', '%+02.6f', '%6i', '%6i', '%6i', '%6i', '%6i', '%7s']

def night_from_filename(filename):
    filename = os.path.basename(filename)
    night = os.path.splitext(filename)[0]

    return night

def compare_fred_data(A, B):
    """Compares two FRED entries stored as numpy structured arrays.
    Returns True if they match identically, False otherwise."""

    same_point = True

    for key in fred_col_names:
        if A[key] != B[key]:
            same_point = False

    return same_point

def read_fred_manual(filename):
    """Reads and tokenizes a FRED file using (simplistic) parsing methods.
    Lines with errors will generate a warning message.

    This function handles errors in the FRED files more gracefully than
    read_fred, but is significantly slower.
    """

    # define a flags
    flag_missing_columns = False

    data = []
    num_fred_cols = len(fred_col_names)
    dtype={'names': fred_col_names,'formats': fred_col_types}

    with open(filename, 'r') as infile:
        for line in infile:
            # remove whitespace from the lines
            line = line.strip()

            # skip comment lines
            if line[0] == "#":
                continue

            # tokenize the line
            line = line.split()

            if len(line) < num_fred_cols:
                flag_missing_columns = True
                continue

            # initial checks are ok. Append the data.
            data.append(line)

    if len(data) == 0:
        return None

    if flag_missing_columns:
        print("WARNING: %s is missing required columns" % (filename))

    # Convert the data into a numpy array:
    data = np.asarray(data, dtype=dtype)
    return data

def read_fred(filename):
    """Reads in an APASS FRED file and returns the result as a numpy structured
    array with columns as specified in fred_col_names plus fredbin.fredbin_extra_cols.
    In rare situations, this function can return None."""
    dtype={'names': fred_col_names,'formats': fred_col_types}
    data = []

    # first attempt to load using numpy.loadtext
    try:
        data = np.loadtxt(filename, dtype=dtype)
    except (ValueError,IndexError):
        # if this fails, try parsing manually
        data = read_fred_manual(filename)

    num_rows = len(data)

    night_name = night_from_filename(filename)
    night_names = [night_name] * num_rows

    # append extra type columns for zone_id, node_id, and container_id
    tmp = np.zeros(num_rows)
    data = nprf.append_fields(data, fredbin_extra_col_names, [tmp, tmp, tmp, night_names],
                              dtypes=fredbin_extra_col_types)

    return data


def read_fredbin(filename):
    """Reads in an APASS .fredbin file."""

    dtype={'names': fredbin_col_names,'formats': fredbin_col_types}
    data = np.fromfile(filename, dtype)

    return data

def write_txt(fredbin_data, filename):
    header  = "APASS .fredbin text file dump. Format is as follows:\n"
    header += "Column names: "   + ','.join(fredbin_col_names) + "\n"
    header += "Column types: "   + ','.join(fredbin_col_types) + "\n"
    header += "Column formats: " + ','.join(fredbin_col_fmt)   + "\n"
    np.savetxt(filename + ".txt", fredbin_data, header=header, fmt=fredbin_col_fmt)
