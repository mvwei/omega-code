import numpy as np
import pandas
import matplotlib.pyplot as plt

from functions import great_circle_distance
from classes import CoordinateField

def get_lat_degree_interval(dy, r=6367.0):
    """
    We're doing arcs along a circle here. Woo!

    arclength = r * theta, where theta is angle of difference in radians.
    """
    return np.degrees(dy/r)

def get_lon_degree_interval(lat, dx, r=6367.0):
    """
    No change in latitude. All we need to do is figure the length of
    the "radius" of the slice of the sphere at the provided latitude.

    rlat = r * cos(lat)
    dx = rlat * theta = r * cos(lat) * theta
    """
    rad_lat = np.radians(lat)

    return np.degrees(dx/(r * np.cos(rad_lat)))


def get_equidistant_lat_lon(north_lat_bound, south_lat_bound, west_lon_bound, east_lon_bound, dx, dy):
    r = 6367

    lat_interval = get_lat_degree_interval(dy)
    lat_array = np.arange(south_lat_bound, north_lat_bound, lat_interval)

    dlon0 = get_lon_degree_interval(lat=south_lat_bound, dx=dx)
    nlon0 = int(np.ceil((east_lon_bound - west_lon_bound) / dlon0))

    lon_array = [
        (get_lon_degree_interval(lat, dx) * np.arange(0, nlon0) + west_lon_bound)
        for lat in lat_array
    ]

    lat_array_2d = np.tile(lat_array, (nlon0, 1)).transpose()

    return lat_array_2d, lon_array

def visualize(data):
    """
    I don't know, what the fuck do I know?
    """
    return

def get_2d_matrix_row(i, j, k, nx, ny):
    return k*nx*ny + i*nx + j
