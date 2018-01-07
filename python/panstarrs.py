# Functions for reading Pan-STARRS database queries.

# system includes
import numpy as np
import numpy.lib.recfunctions as nprf

def read_panstarrs(filename):
    """Reads Pan-STARRS data from CSV database query dumps

    Returns it as a numpy structured array with columns as specified below.
    """


    col_names = ['objID', 'nDetections', 'ra', 'dec',
                 'sg', 'sg_sig', 'sr', 'sr_sig', 'si', 'si_sig', 'sz', 'sz_sig']
    col_types = ['int64', 'int', 'float32', 'float32',
                 'float32', 'float32', 'float32', 'float32',
                 'float32', 'float32','float32', 'float32']

    dtype={'names': col_names,'formats': col_types}
    data = np.loadtxt(filename, skiprows=2, delimiter=',', dtype=dtype)

    # append extra type columns for zone_id, node_id, and container_id
    extra_col_names = ['zone_id', 'node_id', 'container_id']
    extra_col_types = ['int', 'int', 'int']
    tmp = np.zeros(len(data))
    data = nprf.append_fields(data, extra_col_names, [tmp, tmp, tmp],
                              dtypes=extra_col_types)

    # filter out stars with large uncertainties
    indexes = np.where((data['sg_sig'] < 999) &
                       (data['sr_sig'] < 999) &
                       (data['si_sig'] < 999) &
                       (data['sz_sig'] < 999))

    data = data[indexes]

    return data
