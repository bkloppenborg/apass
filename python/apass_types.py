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
        self.num_data = 0

        dr = 1. / (60 * 60) # 1 arcsecond in degrees
        dx = dr / cos(y * pi /  180)
        dy = dr
        self.rect = Rect(x - dx, x + dx, y - dy, y + dy)
        # self.data is a list of numpy.ndarray objects
        self.data = []
        if data is not None:
            self.append_data(data)

    def __repr__(self):
        """Return a string representation of this object"""
        str_repr = "Zone: %i Node: %i ID: %i, NumData: %i Bounds: [%f %f %f %f]" % (
            self.zone_id, self.node_id, self.container_id, num_data,
            self.rect.x_min, self.rect.x_max, self.rect.y_min, self.rect.y_max)

        return str_repr

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

        container.num_data = dict_['num_data']
        container.rect = Rect.from_dict(dict_['rect'])
        return container

    def append_data(self, data):
        """Appends the specified data to this object"""
        self.data.append(data)
        self.num_data = len(self.data)

    def merge(self, other):
        """Merges two RectContainer Instances, growing their bounding rectangles
        other -- another RectContainer instance"""

        # conditionally reflect the data about RA 0 <-> 360
        rect, data = apass.reflect_rect_and_data(other.rect, other.data)

        # grow the rectangle
        self.rect.expand(rect)

        # re-number the other container's data IDs and append it to this node's data
        for i in range(0, len(data)):
            data[i]['node_id'] = self.node_id
            data[i]['container_id'] = self.container_id

        self.data.extend(data)
        other.data = []

        self.num_data = len(self.data)

    def overlaps(self, other):
        """Determines if this RectContainer overlaps with another RectContainer

        other -- another RectContainer instance"""
        overlaps = False

        # If both nodes span the north or south pole, they automatically overlap
        if self.rect.y_max > 90 and other.rect.y_max > 90:
            overlaps = True
        elif self.rect.y_min < -90 and other.rect.y_min < -90:
            overlaps = True
        else:
            s_rect = apass.reflect_rect(self.rect)
            o_rect = apass.reflect_rect(other.rect)
            overlaps = s_rect.overlaps(o_rect)

        return overlaps

    def get_data(self):
        """Returns this container's data as a list."""

        # ensure the data is marked as belonging to this container.
        for datum in self.data:
            datum['zone_id'] = self.zone_id
            datum['node_id'] = self.node_id
            datum['container_id'] = self.container_id

        return self.data

    def clear_data(self):
        """Erases all data contained in this node"""
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

    def get_overlapping_containers(self, other_container):
        """Returns a list of containers that overlaps with other_container.
        If remove=True (default), any overlapping containers will be removed from
        this container, otherwise they will be left intact."""
        containers = []
        for container in self.containers:
            if container.overlaps(other_container):
                containers.append(container)

        return containers

    def remove_container(self, container):
        self.containers.remove(container)

    def remove_containers(self, containers):
        for container in containers:
            self.remove_container(container)

    def insert(self, x, y, data):
        """Stores the data inside of a container encapsulated by this node."""

        # create a container to store this data
        container = RectContainer(x, y, data)

        # find any overlapping containers, sort by number of contained data in
        # descending order
        overlappers = self.get_overlapping_containers(container)
        overlappers.sort(key=lambda x: x.num_data, reverse=True)

        # Either append the container, or merge things together
        if len(overlappers) == 0:
            self.containers.append(container)
        else:
            bigger_cont = overlappers[0]
            other_conts = overlappers[1:]

            bigger_cont.merge(container)
            for other in other_conts:
                bigger_cont.merge(other)

            self.remove_containers(other_conts)

    def insert_or_drop(self, x, y, data, distance=1):
        """Stores data inside of a container encapsulated by this node if a suitable
        container is found within the specified distance"""

        container = RectContainer(x, y, data)

        # find any overlapping containers, sort by number of contained data in
        # descending order
        overlappers = self.get_overlapping_containers(container)
        overlappers.sort(key=lambda x: x.num_data, reverse=True)

        # Either append the container, or merge things together
        if len(overlappers) == 0:
            pass # drop the entry
        else:
            bigger_cont = overlappers[0]
            other_conts = overlappers[1:]

            bigger_cont.merge(container)
            for other in other_conts:
                bigger_cont.merge(other)

            self.remove_containers(other_conts)

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

        for i in range(0, len(self.containers)):
            container = self.containers[i]
            # update the identifiers
            container.zone_id = self.zone_id
            container.node_id = self.node_id
            container.container_id = i

            container.save_data(filehandle)

