import numpy as np

from functions import (
    first_derivative,
    second_derivative,
)

class ScalarField:
    """
    Given a N dimensional numpy array, here's a class that will
    ideally give us everything we could ever want. For derivatives.

    DIMENSIONS: could be (time, plvl, lat, lon), (plvl, lat, lon), OR
    (lat, lon). Who knows.
    """
    def __init__(self, values, dx, dy):
        """

        """
        self.values = values

        self.dx = dx
        self.dy = dy

        self.ddx = None
        self.ddy = None

        self.laplacian = None

    def get_ddx(self):
        if self.ddx is not None:
            return self.ddx

        self.ddx = first_derivative(f=self.values, axis=-1, distance=self.dx)

        return self.ddx

    def get_ddy(self):
        if self.ddy is not None:
            return self.ddy

        self.ddy = first_derivative(f=self.values, axis=-2, distance=self.dy)

        return self.ddy

    def get_laplacian(self):
        if self.laplacian is not None:
            return self.laplacian

        d2dx2 = second_derivative(f=self.values, axis=-1, distance=self.dx)
        d2dy2 = second_derivative(f=self.values, axis=-2, distance=self.dy)

        self.laplacian = self.d2dx2 + self.d2dy2

        return self.laplacian

class VectorField:
    """
    Right now, let's just assume we're using a two-dimensional vector
    field. No need to get pointlessly fancy.
    """
    def __init__(self, x_values, y_values, dx, dy):
        self.x_values = x_values
        self.y_values = y_values

        self.dx = dx
        self.dy = dy

    def divergence(self):
        dXdx = first_derivative(f=self.x_values, axis=-1, distance=self.dx)
        dYdy = first_derivative(f=self.y_values, axis=-2, distance=self.dy)

        return dXdx + dYdy


class CoordinateField:
    """
    REMINDER: y is dimension 0, x is dimension 1.
    """
    def __init__(self, lat, lon, distance_func, levels):
        """
        lat: 1D array of latitudes.
        lon: 1d array of longitudes.
        lat x lon makes up our entire grid.
        """
        self.nx = lon.size
        self.ny = lat.size

        self.distance_func = distance_func

        self.levels = levels

        lon2d = np.tile(lon, (self.ny, 1))
        lat2d = np.tile(lat, (self.nx, 1)).transpose()

        self.coors = np.dstack((lat2d, lon2d))

        self.dx = self._get_dx()
        self.dy = self._get_dy()
        self.dp = self._get_dp()

    def _get_dx(self):
        """
        Returns a 1D array of distances (km) between horizontal points.
        For a given latitude, the distance between points should be the same.
        Horizontal distance decreases as latitude increases, because globes
        and stuff.
        """
        slice1 = self.coors[:, :-1]
        slice2 = self.coors[:, 1:]

        distances = self.distance_func(slice1, slice2)

        distance_arr = distances[:, 0]

        # if all is well, then the distances between points should be the same
        # along constant latitude.
        # Minus one because we're calculating the distance between points, and
        # therefore we have one fewer distances than we have points.
        expected_distances = np.tile(distance_arr, (self.nx-1, 1)).transpose()

        if not np.allclose(distances, expected_distances):
            print("Warning: longitudinal distances are not constant along constant latitude. dx is not set.")
            return

        return distance_arr

    def _get_dy(self):
        """
        Returns a constant representing the latitudinal distance between points of
        constant longitude. Result should be constant for all points, because
        we're increasing a constant number of degrees per row.
        """
        slice1 = self.coors[:-1, :]
        slice2 = self.coors[1:, :]

        distances = self.distance_func(slice1, slice2)

        expected_distance = distances[0, 0]

        if not np.allclose(distances, expected_distance):
            print("Warning: latitudinal distances are not constant along constant longitude. dy is not set.")
            return

        return expected_distance

    def _get_dp(self):
        """
        This is mostly a sanity check to see that all the pressure levels
        are evently distributed.
        """
        first = self.levels[0]
        last = self.levels[-1]

        evenly_spaced = np.linspace(first, last, num=self.levels.size)

        if not np.array_equal(self.levels, evenly_spaced):
            print("Warning: pressure levels are not evenly spaced. Please fix.")
            return

        return self.levels[1] - self.levels[0]
