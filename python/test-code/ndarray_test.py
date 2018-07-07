#!/usr/bin/python

import numpy as np


fred_col_names  = ['ra',       'dec']
fred_col_types  = ['float64',  'float64']
dtype={'names': fred_col_names,'formats': fred_col_types}

print("generating and saving data")
filename = '/tmp/test'
outfile = open(filename, 'wb')
for i in range(0, 10):
    values = (1.0 * i, 2.0 * i)
    data = np.array(values, dtype=dtype)
    print data
    np.save(outfile, data, allow_pickle=False)

    print i, values

outfile.close()

print("Reading data from file")
infile = open(filename, 'r')
data = None
while(True):

    try:
        temp = np.load(infile)
    except IOError:
        break

    print temp

    if data is None:
        data = temp
    else:
        data = np.append(data, temp)

print("Final data array")
print data
print data['ra']
print data.dtype
