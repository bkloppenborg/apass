from quadtree import *
from quadtree_types import *
from math import cos, pi
import apass
import os

class RectContainer(dict):
    """A class which implements a storage container with rectangular boundaries."""

    def __init__(self, x, y, data, zone_id=-1, node_id=-1, container_id=-1):
        dict.__init__(self)
        self.__dict__ = self

        self.zone_id = zone_id
        self.node_id = node_id
        self.container_id = container_id
        self.moved = False
        self.moved_zone_id = -1
        self.moved_node_id = -1
        self.moved_container_id = -1
        self.num_data = 0

        dr = 2. / (60 * 60) # 2 arcsecond in degrees
        dx = dr * cos(y * pi /  180)
        dy = dr
        self.rect = Rect(x - dx, x + dx, y - dy, y + dy)
        # self.data is a list of numpy.ndarray objects
        self.data = []
        if data is not None:
            self.append_data(data)

    def contains(self, x, y):
        """Determines if the point (x,y) is contained within this container's bounds"""
        return self.rect.contains(x,y)

    @staticmethod
    def from_dict(dict_):
        """Restores the object and its metadata.
        To restore the data, call load_from_directory()"""
        zone_id = dict_['zone_id']
        node_id = dict_['node_id']
        container_id = dict_['container_id']

        container = RectContainer(0, 0, None,
                                  zone_id = zone_id, node_id = node_id, container_id = container_id)

        container.moved = dict_['moved']
        container.moved_zone_id = dict_['moved_zone_id']
        container.moved_node_id = dict_['moved_node_id']
        container.moved_container_id = dict_['moved_container_id']
        container.num_data = dict_['num_data']
        container.rect = Rect.from_dict(dict_['rect'])
        return container

    def append_data(self, data):
        """Appends the specified data to this object"""
        self.data.append(data)
        self.num_data += 1

    def merge(self, other, mark_moved=False):
        """Merges two RectContainer Instances, growing their bounding rectangles
        other -- another RectContainer instance"""

        # grow the border:
        self.rect.expand(other.rect)

        # re-number the other container's data IDs and append it to this node's data
        for i in range(0, len(other.data)):
            other.data[i]['node_id'] = self.node_id
            other.data[i]['container_id'] = self.container_id
        self.data.extend(other.data)
        other.data = []

        # mark the other container's data as moved to this container
        if mark_moved:
            other.moved = True
            other.moved_zone_id = self.zone_id
            other.moved_node_id = self.node_id
            other.moved_container_id = self.container_id

    def overlaps(self, other):
        """Determines if this RectContainer overlaps with another RectContainer

        other -- another RectContainer instance"""
        return self.rect.overlaps(other.rect)

    def save_data(self, filehandle):
        """Writes this container's data to the specified filehandle"""

        if len(self.data) == 0:
            return

        for datum in self.data:
            filehandle.write(datum)
        self.data = []

    def get_corners(self):
        return self.rect.get_corners()

class RectLeaf(QuadTreeNode):
    """A class for a QuadTree leaf that contains a list of rectangles with
    data."""

    def __init__(self, rect, depth, parent=None, zone_id=-1, node_id=-1):
        QuadTreeNode.__init__(self, rect, depth, parent)

        self.zone_id = zone_id
        self.node_id = node_id
        self.containers = []

    @staticmethod
    def from_dict(rect, depth, dict_):
        """Restores a RectLeaf. See QuadTreeNode.from_dict"""
        zone_id = dict_['zone_id']
        node_id = dict_['node_id']

        # restore the RectLeaf node
        leaf = RectLeaf(rect, depth, zone_id=zone_id, node_id=node_id)
        # restore any RectContainers (metadata only)
        for container_info in dict_['containers']:
            leaf.containers.append(RectContainer.from_dict(container_info))

        return leaf

    def get_container(self, x, y):
        """Returns a reference to the container which contains (x,y) or None if
        no such container exists"""
        for container in self.containers:
            if container.contains(x,y):
                return container

        return None

    def get_overlapping_containers(self, other_container, remove=False):
        """Returns a list of containers that overlaps with other_container.
        If remove=True (default), any overlapping containers will be removed from
        this container, otherwise they will be left intact."""
        containers = []
        for container in self.containers:
            if container.overlaps(other_container):
                containers.append(container)

        # remove the containers from self.containers if requested
        if remove:
            for container in containers:
                self.containers.remove(container)

        return containers

    def insert(self, x, y, data):
        """Stores the data inside of a container encapsulated by this node."""

        # build a rectangle container for this data
        container = RectContainer(x, y, data)
        overlappers = self.get_overlapping_containers(container, remove=True)
        for other in overlappers:
            container.merge(other)
        self.containers.append(container)

    def load_data(self, data):
        """Restores the specified data to the container. Used in object restoration."""
        if self.zone_id < 0 and self.node_id < 0:
            raise RuntimeError("Cannot load data for an uninitialized node!")

        # build a dictionary of container ID -> containers
        container_dict = dict()
        for container in self.containers:
            container_id = container.container_id
            container_dict[container_id] = container

        # insert the  data
        for i in range(0, len(data)):
            container_id = data[i]['container_id']
            container_dict[container_id].append_data(data[i])

    def save_data(self, filehandle):
        """Saves the data in this node to the file handle savefile"""
        for container in self.containers:
            container.save_data(filehandle)
