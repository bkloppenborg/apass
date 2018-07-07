#!/bin/python

# system includes
import sys, os
import timeit
import numpy as np

import sys, os
sys.path.append(os.path.join(sys.path[0],'../'))

import fred

def read_fredbin_A(filename):
    # Reads data as numpy structured array, converts to list, concatenates
    # and returns numpy structured array.

    infile = open(filename, 'rb')
    data = []
    while(True):

        try:
            temp = np.load(infile)
            temp = temp.tolist()
            data.extend(temp)
        except IOError:
            break

    # convert the data to a numpy structured array
    data = fred.to_fredbin(data)

def read_fredbin_B(filename):

    # Reads numpy structured array, appends, returns structured array.

    infile = open(filename, 'rb')
    data = None
    while(True):

        try:
            temp = np.load(infile)
        except IOError:
            break

        if data is None:
            data = temp
        else:
            data = np.append(data, temp)

    return data


filename = '/home/data/apass-test/z09199.fredbin'
num_iterations = 10


start = timeit.default_timer()
for i in range(0, num_iterations):
    read_fredbin_A(filename)
stop = timeit.default_timer()
ave_time = (stop-start) / num_iterations

print "read_fred_A time:" + str(ave_time)


start = timeit.default_timer()
for i in range(0, num_iterations):
    read_fredbin_B(filename)
stop = timeit.default_timer()
ave_time = (stop-start) / num_iterations
print "read_fred_B time:" + str(ave_time)
