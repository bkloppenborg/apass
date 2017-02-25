
from functools import partial
import json

class Rect(dict):
    """A class which represents an axis-aligned rectangle"""
    x_min = 0
    x_max = 0
    y_min = 0
    y_max = 0

    def __init__(self, x_min, x_max, y_min, y_max):
        dict.__init__(self)
        self.__dict__ = self
        self.x_min = float(x_min)
        self.x_max = float(x_max)
        self.y_min = float(y_min)
        self.y_max = float(y_max)

    def __repr__(self):
        """Return a string representation of the object"""
        return "[x_min: %f x_max: %f y_min: %f y_max: %f]" % \
            (self.x_min, self.x_max, self.y_min, self.y_max)

    def contains(self, x,y):
        """Determines if the point (x,y) resides within this rectangle.
        The point must be (x,y) >= (x_min,y_min) and (x,y) < (x_max,y_max)"""
        if x >= self.x_min and x < self.x_max and y >= self.y_min and y < self.y_max:
            return True
        return False

    def expand(self, other):
        """Expands this rectangle to include the other rectangle's bounds

        other -- another Rect instance"""
        self.x_min = min(self.x_min, other.x_min)
        self.x_max = max(self.x_max, other.x_max)
        self.y_min = min(self.y_min, other.y_min)
        self.y_max = max(self.y_max, other.y_max)

    @staticmethod
    def from_dict(dict_):
        """Restores a Rect from a dictionary"""
        x_min = dict_['x_min']
        x_max = dict_['x_max']
        y_min = dict_['y_min']
        y_max = dict_['y_max']
        return Rect(x_min, x_max, y_min, y_max)

    def get_corners(self):
        corners = []
        corners.append((self.x_min, self.y_min))
        corners.append((self.x_min, self.y_max))
        corners.append((self.x_max, self.y_min))
        corners.append((self.x_max, self.y_max))
        return corners

    def grow(self, radius):
        """Expands the rectangle by radius in all directions"""
        self.x_min -= radius
        self.x_max += radius
        self.y_min -= radius
        self.y_max += radius

    def overlaps(self, other):
        """Determines if two rectangles overlap"""
        # implemented as "when do two rectangles NOT overlap?"
        if (other.y_max < self.y_min or
            other.x_max < self.x_min or
            other.x_min > self.x_max or
            other.y_min > self.y_max):
            return False
        return True

    def splitIntoQuads(self):
        """Splits a rectangle into four quadrants"""
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

class QuadTreeNode(dict):

    def __init__(self, rect, depth, parent=None):
        """A class representing a quadtree node.

        rect -- the rectangle which this node contains
        depth -- the depth of this node
        parent -- this node's parent"""
        # NOTE: If you add members here, update from_dict's num_keys value!
        dict.__init__(self)
        self.__dict__ = self
        self.parent = parent
        self.children = []
        self.rect = rect
        self.depth = depth

    def __repr__(self):
        """Return a string representation of this object"""
        return "rect: %s, depth: %i" % (str(self.rect), self.depth)

    def contains(self, x, y):
        """Determine if this node contains the point x,y"""
        return self.rect.contains(x,y)

    @staticmethod
    def from_dict(dict_, leafClass=None):
        """Restores a QuadTree node or subclass thereof

        dict -- a dictionary containing this quadtree node's data
        leafClass -- The class which will represent the leaves, must be a
                     subclass of QuadTree"""
        # extract common data for every node
        rect = Rect.from_dict(dict_['rect'])
        depth = dict_['depth']
        json_children = dict_['children']
        parent = dict_['parent'] # always is none

        node = None
        num_keys = len(dict_.keys())
        if num_keys == 4:
            # standard QuadTree node
            node = QuadTreeNode(rect, depth)
            mapfunc = partial(QuadTreeNode.from_dict, leafClass=leafClass)
            node.children = list(map(mapfunc, json_children))
        else:
            # specialized node, use leafClass's from_dict function
            node = leafClass.from_dict(rect, depth, dict_)

        return node

    @staticmethod
    def to_json(tree):
        """Serializes the quadtree tree to a JSON string"""

        # Remove parents to avoid cyclic dependencies
        tree.runFunc(remove_parent_from_children)
        return json.dumps(tree, sort_keys=True, indent=2)

    @staticmethod
    def to_file(tree, filename):
        """Serializes a quadtree to a JSON file"""
        json_str = QuadTreeNode.to_json(tree)
        with open(filename, 'w+') as outfile:
            outfile.write(json_str)

    @staticmethod
    def from_file(filename, leafClass=None):
        """Restores a quadtree from a JSON File

        filename -- the input JSON file
        leafClass -- The class which will represent the leaves, must be a
                     subclass of QuadTree"""
        json_str = open(filename).read()
        return QuadTreeNode.from_json(json_str, leafClass=leafClass)

    @staticmethod
    def from_json(json_str, leafClass=None):
        """Restores the quadtree from a JSON string."""

        # restore the tree from a JSON file
        json_data = json.loads(json_str)
        tree = QuadTreeNode.from_dict(json_data, leafClass=leafClass) # reconstitute
        tree.runFunc(restore_parent_to_children)
        return tree

    def find_leaf(self, x, y):
        """Finds the leaf which contains (x,y). If no such node can be located,
        None is returned."""
        if self.is_leaf():
            return self

        for child in self.children:
            if child.contains(x,y):
                return child.find_leaf(x, y)

        return None

    def find_root(self):
        """Traverse the tree from the leaves to the root and returns the root
        when found. If the tree is broken (i.e. parents are not restored during
        an import operation, the node will return a reference to itself)"""
        if self.parent is not None:
            return parent.find_root()

        return self

    def find_node_containing(self, x, y):
        """Traverse the tree in either direction to find a node which contains
        the pair (x,y). Returns none of no such node can be found."""

        if self.parent is None:
            return None

        if self.parent.contains(x,y):
            return self.parent.find_leaf(x,y)
        else:
            return self.parent.find_node_containing(x,y)

    def get_leaves(self):
        leaves = []
        func = partial(get_leaves, leaves=leaves)
        self.runFunc(func)
        return leaves


    def has_children(self):
        """Determines if this node has children"""
        return len(self.children) > 0

    def insert(self, x, y, datum):
        """Traverse the tree from the root towards the leaves and insert the
        data into the appropriate leaf"""
        for child in self.children:
            if child.contains(x, y):
                return child.insert(x, y, datum)

        # hopefully we never encounter this...
        raise ValueError("Could not find a node containing the point (%f, %f)" % (x,y))

    def is_leaf(self):
        return not self.has_children()

    def is_root(self):
        return self.parent == None

    def runFunc(self, f):
        """Runs the function on the node and its children.
        The function must take a node as a parameter"""
        f(self)

        for child in self.children:
            child.runFunc(f)

    def split_until(self, depth, leafClass=None):
        """Subdivides a quadtree until the specified depth. Leaf nodes
        will be of the  type QuadTreeNode unless leafClass is specified."""

        if(self.depth == depth - 1):
            # split into leaf nodes
            self.split(leafClass=leafClass)

        elif(self.depth < depth):
            # split into normal nodes
            self.split(leafClass=None)

            for child in self.children:
                child.split_until(depth, leafClass=leafClass)

    def size(self):
        """Determine the total number of nodes in the tree"""
        size = 1
        for child in self.children:
            size += child.size()

        return size

    def split(self, leafClass=None):
        """Subdivides a node into four approximately equal sized quadrants"""
        rects = self.rect.splitIntoQuads()
        for rect in rects:
            if leafClass == None:
                self.children.append(QuadTreeNode(rect, self.depth + 1, parent=self))
            else:
                self.children.append(leafClass(rect, self.depth + 1, parent=self))


def remove_parent_from_children(node):
    """Removes this node from being the parent of its children. Used in serialization."""
    for child in node.children:
        child.parent = None

def restore_parent_to_children(node):
    """Restores this node as the parent of its children. Used in de-serialization."""
    for child in node.children:
        child.parent = node

def get_leaves(node, leaves=None):
    if leaves is None:
        raise ValueError("You must specify leaves via keyword using functools.partial")
    if node.is_leaf():
        leaves.append(node)
