
from quadtree import *

class IDLeaf(QuadTreeNode):
    """Class for a QuadTree leaf that points to a file"""

    node_id_counter = 0

    def __init__(self, rect, depth, node_id=None, parent=None):
        QuadTreeNode.__init__(self, rect, depth, parent)

        self.node_id = node_id
        if self.node_id is None:
            self.node_id = IDLeaf.node_id_counter
            IDLeaf.node_id_counter += 1

    @staticmethod
    def from_dict(rect, depth, dict_):
        """Restores a IDLeaf. See QuadTreeNode.from_dict"""
        node_id = dict_['node_id']
        return IDLeaf(rect, depth, node_id=node_id)

    def insert(self, x, y, data):
        """Returns an ID into which insert should occur.
        NOTE: DOES NOT INSERT ANY DATA"""
        return self.node_id


