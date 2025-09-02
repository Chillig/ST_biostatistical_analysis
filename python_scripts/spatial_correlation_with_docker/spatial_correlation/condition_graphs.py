"""Define nearest neighbor spots and connect them via a graph
    File name: condition_graphs.py
    Author: Christina Hillig
    Date created: 2/01/2021
    Date last modified: 09/01/2025
    Python Version: 3.8
"""


# Graphs
import networkx
from networkx.algorithms.components.connected import connected_components

# For cKDTree / KDTree
from sklearn.neighbors import KDTree
from scipy import spatial

import numpy as np


class SpotConnectivity:
    """
    Class to build connectivity graphs and clusters of spatial points
    using KDTree or cKDTree.

    Methods
    -------
    to_graph(points)
        Create a graph from nodes (spots) and edges.
    to_edges(points)
        Convert a sequence of points into edge pairs.
    kdtree_clustering(coordinates, distance)
        Find nearest neighbor points using KDTree within a distance.
    get_connectedspots(coordinates, distance)
        Use cKDTree to find conditionally connected spot pairs.
    get_connections(coordinates, distance, tree_type)
        Wrapper to return connected components or spot pairs depending on tree_type.
    """

    def __init__(self, coordinates: np.ndarray):
        """
        Parameters
        ----------
        coordinates : np.ndarray
            Array of shape (n_samples, 2) with x, y positions of spots.
        """
        self.coordinates = np.array(coordinates)

    @staticmethod
    def to_graph(points):
        """Create a graph from nodes and edges inferred from sequential points."""
        graph = networkx.Graph()
        for part in points:
            graph.add_nodes_from(part)
            graph.add_edges_from(SpotConnectivity.to_edges(part))
        return graph

    @staticmethod
    def to_edges(points):
        """Yield consecutive edges from a list of points."""
        it = iter(points)
        last = next(it)
        for current in it:
            yield last, current
            last = current

    def kdtree_clustering(self, distance: float = 2.0):
        """
        Find nearest neighbor points using KDTree within max distance.
        """
        kdt = KDTree(self.coordinates, leaf_size=30, metric='euclidean')
        dist_points, ind_points = kdt.query(self.coordinates, k=len(self.coordinates), return_distance=True)

        m_distance = dist_points <= distance
        nn_points = []
        for ind, combo in enumerate(ind_points):
            temp_nnpoints = list(combo[m_distance[ind]])
            nn_points.append(temp_nnpoints)

        return nn_points

    def get_connectedspots(self, distance: float):
        """
        Find connected spot pairs using cKDTree within a max distance.
        """
        tree = spatial.cKDTree(self.coordinates)
        g = {i: set(tree.query(self.coordinates[i, :], len(self.coordinates))[-1][1:])
             for i in range(self.coordinates.shape[0])}

        connected_spots = []
        for node, candidates in g.items():
            for node2 in candidates:
                if node2 < node:
                    continue
                if node in g[node2] and spatial.distance.euclidean(self.coordinates[node, :],
                                                                   self.coordinates[node2, :]) <= distance:
                    connected_spots.append((node, node2))
        return connected_spots

    def get_connections(self, distance: float = 2.0, tree_type: str = 'KDTree'):
        """
        Get connected components of spots using KDTree or cKDTree.

        Parameters
        ----------
        distance : float
            Maximum distance between spots to be considered connected.
        tree_type : str
            Either 'KDTree' (clustered components) or 'cKDTree' (pairs of connected spots).

        Returns
        -------
        list
            If KDTree: list of connected components (sets of indices).
            If cKDTree: list of connected pairs and their coordinates.
        """
        if tree_type == 'KDTree':
            nn_points = self.kdtree_clustering(distance=distance)
            graph_tree = self.to_graph(nn_points)
            return sorted(connected_components(graph_tree), key=len, reverse=True)
        else:
            connected_spots = self.get_connectedspots(distance=distance)
            lines = [[self.coordinates[i], self.coordinates[j]] for i, j in connected_spots]
            return connected_spots, lines
