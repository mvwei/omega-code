import numpy as np

from functions import (
    get_shifted_by_one,
    horizontal_gradient,
    vertical_gradient,
    laplacian,
    great_circle_distance,
)

# Given a N dimensional numpy array, here's a class that will
# ideally give us everything we could ever want.
class ScalarField:
    def __init__(self, np_array):
        self.values = values

    def ddx(self):
        return

    def ddy(self):
        return

    def ddp(self):
        return

    def d2dx2(self):
        return

    def d2dy2(self):
        return

    def horizontal_gradient(self):
        return

    def laplacian(self):
        return

# Values should be of form (a, b), corresponding to the vector ai + bj.
class VectorField:
    def __init__(self, values):
        self.values = values

    def divergence(self):
        return

# REMINDER: y is dimension 0, x is dimension 1.
# lat: 1D array of latitudes.
# lon: 1d array of longitudes.
# lat x lon makes up our entire grid.
class CoordinateField:
    def __init__(self, lat, lon):
        self.nx = lon.size
        self.ny = lat.size

        self.lon = lon
        self.lat = lat

        lon2d = np.tile(lon, (self.ny, 1))
        lat2d = np.tile(lat, (self.nx, 1)).transpose()

        self.coors = np.dstack((lat2d, lon2d))

        self.dx = None
        self.dy = None

    def get_dx(self):
        """
        Returns a 1D array of distances (km) between horizontal points.
        For a given latitude, the distance between points should be the same.
        Horizontal distance decreases as latitude increases, because globes
        and stuff.
        """
        if self.dx:
            return self.dx

        slice1 = self.coors[:, :-1]
        slice2 = self.coors[:, 1:]

        distances = great_circle_distance(slice1, slice2)

        distance_arr = distances[:, 0]

        # if all is well, then the distances between points should be the same
        # along constant latitude.
        # Minus one because we're calculating the distance between points, and
        # therefore we have one fewer distances than we have points.
        expected_distances = np.tile(distance_arr, (self.nx-1, 1)).transpose()

        if not np.allclose(distances, expected_distances):
            print("Warning: longitudinal distances are not constant along constant latitude. dx is not set.")
            return

        self.dx = distance_arr

        return self.dx

    def get_dy(self):
        """
        Returns a constant representing the latitudinal distance between points of
        constant longitude. Result should be constant for all points, because
        we're increasing a constant number of degrees per row.
        """
        if self.dy:
            return self.dy

        slice1 = self.coors[:-1, :]
        slice2 = self.coors[1:, :]

        distances = great_circle_distance(slice1, slice2)

        expected_distance = distances[0, 0]

        if not np.allclose(distances, expected_distance):
            print("Warning: latitudinal distances are not constant along constant longitude. dy is not set.")
            return

        self.dy = expected_distance

        return self.dy

