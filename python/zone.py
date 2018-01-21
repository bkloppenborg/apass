# file for loading/saving zone data

# system libraries
import apass
import sys
import os

# custom modules
sys.path.append(os.path.join(sys.path[0],'modules', 'FileLock', 'filelock'))
from filelock import FileLock

# APASS-specific things
import apass
from quadtree import QuadTreeNode
from quadtree_types import *
from apass_types import *
from border_info import load_border_info, save_border_info
from fred import read_fredbin

def load_zone_data(tree, save_dir):
    """Restores zone data from the specified save_dir to the tree."""

    # build a nested dictionary that gives us direct access to the
    # underlying containers. Indexing shall be
    #  container = nodes[node_id][container_id]
    node_dict = dict()
    leaves = tree.get_leaves() # these will be of type RectLeaf
    for leaf in leaves:
        node_id = int(leaf.node_id)

        container_dict = dict()
        for container in leaf.containers:
            container_id = int(container.container_id)

            container_dict[container_id] = container

        node_dict[node_id] = container_dict

    # load the data
    zone_id  = leaves[0].zone_id
    filename = save_dir + '/' + apass.name_zone_container_file(zone_id)
    data     = read_fredbin(filename)

    # insert the data *directly* into the container, bipassing normal
    # restoration methods.
    for i in range(0, len(data)):
        datum        = data[i]
        node_id      = int(datum['node_id'])
        container_id = int(datum['container_id'])
        node_dict[node_id][container_id].append_data(datum)

def save_zone_data(tree, directory):
    """Saves zone data from the tree to the specified directory"""

    leaves = tree.get_leaves()
    zone_id = leaves[0].zone_id
    filename = directory + '/' + apass.name_zone_container_file(zone_id)

    # save the rectangle data to file, overwriting the contents
    with open(filename, 'w+b') as data_file:
        for leaf in leaves:
            leaf.save_data(data_file)

def load_zone(save_dir, zone_id):
    """Loads the tree, data, and border info file for the specified zone.
    Returns this data as a dictionary keyed as follows:
     * json_filename
     * container_filename
     * border_info_filename
     * border_info
     * tree
     * loaded
     * lock (a FileLock instance)"""

    #print(" Loading zone %i" % (zone_id))

    zone_dict = dict()
    zone_tree = None
    zone_border_info = None
    load_succeeded = True
    lock = None

    zone_json_file        = save_dir + apass.name_zone_json_file(zone_id)
    zone_container_file   = save_dir + apass.name_zone_container_file(zone_id)
    zone_border_info_file = save_dir + apass.name_zone_border_file(zone_id)

    lock = FileLock(zone_container_file)
    lock.acquire()

    # verify the necessary files exist, if not bail
    if not os.path.isfile(zone_json_file):
#        print(" WARNING: JSON file missing for zone %i" % (zone_id))
        load_succeeded = False
    if not os.path.isfile(zone_border_info_file):
#        print(" WARNING: Border Info file missing for zone %i" % (zone_id))
        load_succeeded = False
    if not os.path.isfile(zone_container_file):
#        print(" WARNING: Container data file missing for zone %i" % (zone_id))
        load_succeeded = False

    if load_succeeded:
        zone_border_info = load_border_info(zone_border_info_file)
        zone_tree = QuadTreeNode.from_file(zone_json_file, leafClass=RectLeaf)
        load_zone_data(zone_tree, save_dir)

    zone_dict['json_filename']        = zone_json_file
    zone_dict['container_filename']   = zone_container_file
    zone_dict['border_info_filename'] = zone_border_info_file
    zone_dict['tree']                 = zone_tree
    zone_dict['border_info']          = zone_border_info
    zone_dict['loaded']               = load_succeeded
    zone_dict['lock']                 = lock
    return zone_dict

def save_zone(save_dir, zone_dict):
    """Saves the zone data to disk. The zone_dict format matches the
    format found in load_zone above."""

    json_file        = zone_dict['json_filename']
    container_file   = zone_dict['container_filename']
    border_info_file = zone_dict['border_info_filename']
    tree             = zone_dict['tree']
    border_info      = zone_dict['border_info']
    load_succeeded   = zone_dict['loaded']
    lock             = zone_dict['lock']

    zone_id = apass.zone_from_name(container_file)
    #print("Saving zone %i" % (zone_id))

    # save the data to disk
    if load_succeeded:
        save_zone_data(tree, save_dir)
        save_border_info(border_info_file, border_info)
        QuadTreeNode.to_file(tree, json_file)

    # release the lock
    lock.release()
