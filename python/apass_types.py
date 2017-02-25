from quadtree import *
from quadtree_types import *
from math import cos, pi
import apass

class RectContainer():

    def __init__(self, x, y, data, zone_id=-1, node_id=-1, container_id=-1):

        dr = 1. / (60 * 60) # 1 arcsecond in degree
        dx = dr * cos(y * pi /  180)
        dy = dr
        self.rect = Rect(x - dx, x + dx, y - dy, y + dy)
        self.data = [data]
        self.zone_id = zone_id
        self.node_id = node_id
        self.container_id = container_id

        #self.__init__(rect, depth, zone_id=zone_id, node_id=node_id, container_id=container_id)

    def __init__self(rect, depth, zone_id=-1,  node_id=-1, container_id=-1):
        self.rect = rect
        self.data = []
        self.zone_id = zone_id
        self.node_id = node_id
        self.container_id = container_id

    def overlaps(self, other):
        """Determines if this RectContainer overlaps with another RectContainer

        other -- another RectContainer instance"""
        return self.rect.overlaps(other.rect)

    def merge(self, other):
        """Merges two RectContainer Instances, growing their bounding rectangles

        other -- another RectContainer instance"""
        self.rect.expand(other.rect)
        self.data.extend(other.data)

    def save(self, directory):

        filename = directory + "/" + \
                   apass.name_rect(self.zone_id, self.node_id, self.container_id)

        outfile = open(filename, 'a+b')
        for datum in self.data:
            outfile.write(datum)
        outfile.close()
        self.data = []

    @staticmethod
    def from_dict(rect, depth, dict_):
        zone_id = dict['zone_id']
        node_id = dict['node_id']
        container_id = dict['container_id']
        return RectContainer(rect, depth,
                             zone_id=zone_id, node_id=node_id, container_id=container_id)


class RectLeaf(QuadTreeNode):
    """A class for a QuadTree leaf that contains a list of rectangles with
    data."""

    def __init__(self, rect, depth, parent=None, node_id=-1):
        QuadTreeNode.__init__(self, rect, depth, parent)

        self.node_id = node_id
        self.containers = []

    def insert(self, x, y, data):
        """Stores the data inside of a container encapsulated by this node."""

        # build a rectangle container for this data
        rc = RectContainer(x, y, data)
        self.merge_containers(rc)
        self.containers.append(rc)

    def merge_containers(self, rc):
        """Finds containers that overlap with the container rc and merges them
        into rc. As such, the merge containers are removed from this object's
        list of containers."""
        # find, merge, and remove containers that overlap with rc
        containers = []
        for container in self.containers:
            if rc.overlaps(container):
                rc.merge(container)
            else:
                containers.append(container)

        self.containers = containers

    @staticmethod
    def from_dict(rect, depth, dict_):
        """Restores a RectLeaf. See QuadTreeNode.from_dict"""
        node_id = dict_['node_id']
        return RectLeaf(rect, depth, node_id=node_id)

