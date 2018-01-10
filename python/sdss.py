# Functions for reading Sloan Digital Sky Survey CSV data

# system includes
import numpy as np
import numpy.lib.recfunctions as nprf

def read_sdss(filename):
    """Reads SDSS CSV formatted data downloaded from the SDSS SkyServer
    http://skyserver.sdss.org/dr14/en/home.aspx
    rectangular search tool.

    Returns it as a numpy structured array with columns as specified below.
    """

    # format is as follows:
    # objid,run,rerun,camcol,field,obj,type,ra,dec,u,g,r,i,z,Err_u,Err_g,Err_r,Err_i,Err_z
    col_names = ['objid', 'run', 'rerun', 'camcol', 'field', 'obj', 'type',
                 'ra', 'dec',
                 'su', 'sg', 'sr', 'si', 'sz',
                 'su_sig', 'sg_sig', 'sr_sig', 'si_sig', 'sz_sig']
    col_types = ['int64', 'int', 'int', 'int', 'int', 'int', 'int',
                 'float32', 'float32',
                 'float32', 'float32', 'float32', 'float32', 'float32',
                 'float32', 'float32', 'float32', 'float32', 'float32']
    dtype={'names': col_names,'formats': col_types}
    data = np.loadtxt(filename, skiprows=2, delimiter=',', dtype=dtype)

    # append extra type columns for zone_id, node_id, and container_id
    extra_col_names = ['zone_id', 'node_id', 'container_id']
    extra_col_types = ['int', 'int', 'int']
    tmp = np.zeros(len(data))
    data = nprf.append_fields(data, extra_col_names, [tmp, tmp, tmp],
                              dtypes=extra_col_types)

    return data
