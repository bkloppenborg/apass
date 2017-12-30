# Functions to read Carlsberg Meridian Circle Catalog data in text format.

# system includes
import numpy as np

def read_cmc14(filename):
    """Reads a CMC14 file and extracts RA, DEC, sr, sr_sig.
    Returns values as a numpy structured array"""

    # the CMC data is in the following format:
    #RA              DEC             filter  mag     mag_err
    #176.917912      0.309793        r'      14.169  0.028

    data = []
    with open(filename, 'r') as infile:
        for line in infile:
            ra     = float(line[0:12])
            dec    = float(line[12:26])
            sr     = float(line[39:46])
            sr_sig = float(line[46:53])

            data.append((ra, dec, sr, sr_sig, 0, 0, 0))

    col_names = ['ra', 'dec', 'sr', 'sr_sig',
                 'zone_id', 'node_id', 'container_id']
    col_types = ['float32', 'float32', 'float32', 'float32',
                 'int', 'int', 'int']
    dtype={'names': col_names,'formats': col_types}
    data = np.asarray(data, dtype=dtype)

    return data
