
from functools import partial
import json

class Rect(dict):
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
        return "[x_min: %f x_max: %f y_min: %f y_max: %f]" % \
            (self.x_min, self.x_max, self.y_min, self.y_max)

    # this function should be abstracted and implemented elsewhere
    def contains(self, x,y):
        if x >= self.x_min and x < self.x_max and y >= self.y_min and y < self.y_max:
            return True
        return False

    @staticmethod
    def from_dict(dict_):
        x_min = dict_['x_min']
        x_max = dict_['x_max']
        y_min = dict_['y_min']
        y_max = dict_['y_max']
        return Rect(x_min, x_max, y_min, y_max)

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

class QuadTreeNode(dict):

    def __init__(self, rect, depth, parent=None):
        # NOTE: If you add members here, update from_dict's num_keys value!
        dict.__init__(self)
        self.__dict__ = self
        self.parent = parent
        self.children = []
        self.rect = rect
        self.depth = depth

    def __repr__(self):
        return "rect: %s, depth: %i" % (str(self.rect), self.depth)

    def contains(self, x, y):
        if self.rect.contains(x,y):
            return True

        return False

    @staticmethod
    def from_dict(dict_, leafClass=None):
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
        """Serializes the quadtree tree to a JSON file."""

        # Remove parents to avoid cyclic dependencies
        tree.runFunc(remove_parent_from_children)
        return json.dumps(tree, sort_keys=True, indent=2)

    @staticmethod
    def to_file(tree, filename):
        json_str = QuadTreeNode.to_json(tree)
        with open(filename, 'w+') as outfile:
            outfile.write(json_str)

    @staticmethod
    def from_file(filename, leafClass=None):
        json_str = open(filename).read()
        return QuadTreeNode.from_json(json_str, leafClass=leafClass)

    @staticmethod
    def from_json(json_str, leafClass=None):
        """Restores the quadtree from a JSON file."""

        # restore the tree from a JSON file
        json_data = json.loads(json_str)
        tree = QuadTreeNode.from_dict(json_data, leafClass=leafClass) # reconstitute
        tree.runFunc(restore_parent_to_children)
        return tree

    def has_children(self):
        return len(self.children) > 0

    def insert(self, x, y, datum):
        for child in self.children:
            if child.contains(x, y):
                child.insert(x, y, datum)

    def is_leaf(self):
        return not self.has_children()

    def is_root(self):
        return self.parent == None

    def runFunc(self, f):
        """Runs the function on the node and its children"""
        f(self)

        for child in self.children:
            child.runFunc(f)

    def split_until(self, depth, leafClass=None):

        if(self.depth == depth - 1):
            # split into leaf nodes
            self.split(leafClass=leafClass)

        elif(self.depth < depth):
            # split into normal nodes
            self.split(leafClass=None)

            for child in self.children:
                child.split_until(depth, leafClass=leafClass)

    def size(self):
        size = 1
        for child in self.children:
            size += child.size()

        return size

    def split(self, leafClass=None):
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
