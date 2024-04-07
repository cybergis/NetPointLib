from sklearn.neighbors import BallTree
import numpy as np

def build_spatial_index(node_coordinates):
    """
    Builds a spatial index for efficient query of nearest nodes.
    
    Args:
        node_coordinates (numpy.ndarray): Array of node coordinates.

    Returns:
        sklearn.neighbors.BallTree: Spatial index for the node coordinates.
    """
    return BallTree(node_coordinates, metric='haversine')

def query_nearest_edges(events, node_coordinates, spatial_index):
    """
    Efficiently finds the nearest graph edge for each event using a spatial index.
    
    Args:
        events (pandas.DataFrame): DataFrame containing event coordinates.
        node_coordinates (numpy.ndarray): Array of node coordinates.
        spatial_index (sklearn.neighbors.BallTree): Spatial index for the nodes.

    Returns:
        list: Indices of the nearest edge for each event.
    """
    # Convert event coordinates to radians for haversine distance
    events_rad = np.radians(events[['y', 'x']].values)
    # Query the spatial index for nearest node
    distances, indices = spatial_index.query(events_rad, return_distance=True)
    return indices.flatten()