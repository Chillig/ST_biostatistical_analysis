"""Define nearest neighbor spots and connect them via a graph
    File name: condition_graphs.py
    Author: Christina Hillig
    Date created: 2/01/2021
    Date last modified: 3/21/2021
    Python Version: 3.7
"""

# Graphs
import networkx
from networkx.algorithms.components.connected import connected_components

# For cKDTree
from sklearn.neighbors import KDTree
from scipy import spatial

import numpy as np


def to_graph(points):
    """Create a graph from nodes (spots) and edges (connections between spots in a graph)

    Parameters
    ----------
    points : list of tuples

    Returns
    -------

    """
    graph = networkx.Graph()
    for part in points:
        # each sublist is a bunch of nodes
        graph.add_nodes_from(part)
        # it also implies a number of edges:
        graph.add_edges_from(to_edges(part))
    return graph


def to_edges(points):
    """treat 'points' as a Graph and returns it's edges

    Parameters
    ----------
    points

    Returns
    -------
    to_edges(['a','b','c','d']) -> [(a,b), (b,c),(c,d)]

    """
    it = iter(points)
    last = next(it)

    for current in it:
        yield last, current
        last = current


def kdtree_clustering(coordinates, distance=2.0):
    """Use KDTree for fast generalized N-point problems to find nearest neighbor points

    Parameters
    ----------
    coordinates : array of tuples
        [(x1, y2), (x2, y2), ..]
    distance : float
        maximal distance to next point

    Returns
    -------

    """
    # USE KDTree to determine connected cyto+ spots with a distance d
    kdt = KDTree(np.array(coordinates), leaf_size=30, metric='euclidean')
    dist_points, ind_points = kdt.query(np.array(coordinates), k=len(np.array(coordinates)), return_distance=True)

    # Get points which are closely located (<= distance)
    m_distance = dist_points <= distance
    nn_points = []
    for ind, combo in enumerate(ind_points):
        temp_nnpoints = list(combo[m_distance[ind]])
        nn_points.append(temp_nnpoints)

    # return nearest neighbor spots
    return nn_points


def get_connectedspots(coordinates, distance):
    """

    Parameters
    ----------
    coordinates : array of tuples
    distance : float

    Returns
    -------

    """
    tree = spatial.cKDTree(coordinates)
    # each key in g points to indices of 6 nearest cyto+ spots
    g = {i: set(tree.query(coordinates[i, :], len(coordinates))[-1][1:]) for i in range(coordinates.shape[0])}

    connected_spots = list()
    for node, candidates in g.items():
        for node2 in candidates:
            if node2 < node:
                # avoid double-counting
                continue
            if node in g[node2] and spatial.distance.euclidean(coordinates[node, :], coordinates[node2, :]) <= distance:
                connected_spots.append((node, node2))

    return connected_spots


def get_connections(coordinates, distance=2, tree_type='KDTree'):
    """Get spots which are conditionally connected

    Parameters
    ----------
    coordinates : array of tuples
    distance : float
    tree_type : str

    Returns
    -------

    """
    if tree_type == 'KDTree':
        nn_points = kdtree_clustering(coordinates=coordinates, distance=distance)
        graph_tree = to_graph(nn_points)
        return sorted(connected_components(graph_tree), key=len, reverse=True)
    else:
        # USE cKDTree to determine connected cyto+ spots with a distance of d <= 2
        conneted_spots = get_connectedspots(coordinates=np.array(coordinates), distance=1.5)
        lines = [[np.array(coordinates)[i], np.array(coordinates)[j]] for i, j in conneted_spots]
        return conneted_spots, lines
