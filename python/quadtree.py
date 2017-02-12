
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
        return "[x_min: %f x_max: %f y_min: %f y_max: %f]" % \
            (self.x_min, self.x_max, self.y_min, self.y_max)

class QuadTreeNode():

    def __init__(self, rect, depth, parent=None):
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


    def split(self, leafClass=None):
        rects = self.rect.splitIntoQuads()
        for rect in rects:
            if leafClass == None:
                self.children.append(QuadTreeNode(rect, self.depth + 1, parent=self))
            else:
                self.children.append(leafClass(rect, self.depth + 1, parent=self))

