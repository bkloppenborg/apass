
from quadtree import *

class FileStoreLeaf(QuadTreeNode):
    """Class for a QuadTree leaf that points to a file"""

    fileCounter = 0
    filestores = None

    def __init__(self, rect, depth, file_id=None, parent=None, filestores=None):
        QuadTreeNode.__init__(self, rect, depth, parent)
        self.filestores = filestores

        self.file_id = file_id;
        if self.file_id is None:
            self.file_id = fileCounter
            fileCounter += 1

    @staticmethod
    def from_dict(rect, depth, dict_):
        file_id = dict_['file_id']
        return FileStoreLeaf(rect, depth, file_id=file_id)

    def insert(self, x, y, data):
        """Inserts the data into the corresponding file"""
        self.filestores.insert(self.file_id, data)
