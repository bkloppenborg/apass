# Functions to make, read, and write border info files.

# system libraries
import json

# APASS-specific
from apass import name_container

def make_border_info(container):
    """Creates a border info entry for the specified container."""
    info = dict()
    zone_id      = container.zone_id
    node_id      = container.node_id
    container_id = container.container_id

    name = name_container(zone_id, node_id, container_id)

    info['name']         = name
    info['zone_id']      = zone_id
    info['node_id']      = node_id
    info['container_id'] = container_id
    info['center']       = container.rect.get_center()

    nested = {name: info}
    return nested

def save_border_info(filename, information):
    """Saves the border information to the specified file."""
    json_str = json.dumps(information)
    with open(filename, 'w') as outfile:
        outfile.write(json_str)

def load_border_info(filename):
    """Loads the border information from the specified file"""
    json_str = open(filename).read()
    info = json.loads(json_str)
    return info
