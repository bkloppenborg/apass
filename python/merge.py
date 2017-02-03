#!/bin/python

import argparse
import numpy as np

from apass import *
import matplotlib.pyplot as plt
import os.path as path

node_id = 0;

class Rect():
    x_min = 0
    x_max = 0
    y_min = 0
    y_max = 0

    def __init__(self, x_min, x_max, y_min, y_max):
        self.x_min = float(x_min)
        self.x_max = float(x_max)
        self.y_min = float(y_min)
        self.y_max = float(y_max)

    def splitIntoQuads(self):
        output = []
        width =  self.x_max - self.x_min
        height = self.y_max - self.y_min
        dw = width / 2
        dh = height / 2
        x_c = self.x_min + dw
        y_c = self.y_min + dh

        output.append(Rect(self.x_min, x_c, self.y_min, y_c)) # tl
        output.append(Rect(x_c, self.x_max, self.y_min, y_c)) # tr
        output.append(Rect(self.x_min, x_c, y_c, self.y_max)) # bl
        output.append(Rect(x_c, self.x_max, y_c, self.y_max)) # br

        return output

    # this function should be abstracted and implemented elsewhere
    def contains(self, x,y):
        if x >= self.x_min and x < self.x_max and y >= self.y_min and y < self.y_max:
            return True
        return False

    def __repr__(self):
        return "[x_min: %f x_max: %f y_min: %f y_max: %f]" % (self.x_min, self.x_max, self.y_min, self.y_max)

class QuadTree():
    leaves = []
    all_nodes = []
    root_node = None
    bucket_size = 10

    def __init__(self, bounds, bucket_size):
        self.root_node = Node(bounds, 0, 0, bucket_size = bucket_size)

    def insert(self, datum):
        self.root_node.insert(datum)

    def size(self):
        return self.root_node.size()
class NodeData():

    def __init__(self, x, y, index):
        self.x = float(x)
        self.y = float(y)
        self.index = index;

    def __repr__(self):
        return "x: %f y: %f index: %i" % (self.x, self.y, self.index)

class Node():

    def __init__(self, rect, depth, id, bucket_size = 10):
        self.data = []
        self.has_children = False
        self.rect = rect
        self.depth = depth
        self.node_id = id
        self.threshold = bucket_size

    def insert(self, data):
        if self.has_children:
            for child in self.data:
                if(child.contains(data)):
                    child.insert(data)

        else:
            self.data.append(data)
            if len(self.data) + 1 > self.threshold:
                self.split()


    def contains(self, data):
        return self.rect.contains(data.x, data.y)

    def split(self):
        rects = self.rect.splitIntoQuads()
        children = []
        for rect in rects:
            global node_id
            node_id += 1
            children.append(Node(rect, self.depth + 1, node_id))

        for datum in self.data:
            for child in children:
                if child.contains(datum):
                    child.insert(datum)

        self.data = children
        self.has_children = True

    def size(self):
        N = 0

        if self.has_children:
            for child in self.data:
                N += child.size()

        else:
            N = len(self.data)

        return N

    def __repr__(self):
        return "rect: %s, N: " % (str(self.rect), self.data.size)

def open_memmap(filename, dtypes, existing_size, new_size):
    # numpy.memmap containers can't be resized using traditional functions,
    # instead we use this workaround:
    # http://stackoverflow.com/questions/20932361/resizing-numpy-memmap-arrays
    mode = 'r+'
    if existing_size == 0:
        mode = 'w+'

    return np.memmap(filename, dtype=dtypes, mode=mode, shape=(new_size,))


def main():

    parser = argparse.ArgumentParser(description='Merges .fred photometric files')
    parser.add_argument('input', nargs='+', help="Input files which will be merged")

    args = parser.parse_args()

    outdir = '/home/data/apass-source/scratch/'
    output = outdir + 'test.out'
    inputs = args.input
    mmap_file = './temp.bin'

    bounds = Rect(0, 360, -90, 90)
    bucket_size = 10000
    tree = QuadTree(bounds, bucket_size)

    index = 0;
    dtype={'names': apass_col_names,'formats': apass_col_types}

    for filename in inputs:
        print("Reading: " + filename)
        # read in the data file
        data = np.loadtxt(filename, dtype=dtype)

        print("Finished Reading " + filename)

        print("Inserting data into memmap")
        # numpy.memmap objects are a bit buggy in Python, so we open, flush,
        # and close the file every time we append
        #mmap_data = open_memmap(mmap_file, dtype, index, index + data.size)
        #mmap_data[index:] = data
        #mmap_data.flush()
        #del mmap_data

        print("Inserting data into QuadTree")

        for datum in np.nditer(data):
            temp = NodeData(datum['ra'], datum['dec'], index)
            tree.insert(temp)
            index += 1

        #plt.scatter(data['ra'], data['dec'])
        #plt.show()

    print(index)
    print(tree.size())


if __name__ == "__main__":
    main()
