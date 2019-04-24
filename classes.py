import numpy as np

from functions import (
    get_shifted_by_one,
    horizontal_gradient,
    vertical_gradient,
    laplacian,
    great_circle_distance,
)

# Given a three-dimensional NumPy array, make a class around this. This is
# so that we can easily access derivative information, so that we can figure
# out how to write some readable equations. We definitely need some readable
# equations.
class ScalarField:
    def __init__(self, np_array):
        self.values = values

    def ddx(self):
        return

    def ddy(self):
        return

    def ddp(self):
        return

    def gradient(self):
        return

    def laplacian(self):
        return

# Values should be of form (a, b), corresponding to the vector ai + bj
class VectorField:
    def __init__(self, values):
        self.values = values

    def divergence(self):
        return

class CoordinateField:
    def __init__(self, lat, lon):
        self.coors = np.dstack((lat, lon))

        xm1, xp1, ym1, yp1 = get_shifted_by_one(self.coors)

        self.xm1 = xm1
        self.xp1 = xp1
        self.ym1 = ym1
        self.yp1 = yp1

    def get_dx(self):
        """
        Should return a field that represents dx at that specific point
        while calculating a normal x gradient. For interior points, we calculate
        dx as the distance between x+1 and x-1. On the left edge, we calculate
        dx as the difference between x and x+1. On the right edge, we calculate
        dx as x and x-1.
        """
        if hasattr(self, "dx"):
            return self.dx

        self.dx = great_circle_distance(self.xm1, self.xp1)
        return self.dx

    def get_dy(self):
        """
        Should return a field that represents dx at that specific point
        while calculating a normal y gradient. For interior points, we calculate
        dy as the distance between y+1 and y-1. On the "north" edge, we calculate
        dy as the difference between y and y+1. On the  edge, we calculate
        dy as y and y-1.
        """
        if hasattr(self, "dy"):
            return self.dy

        self.dy = great_circle_distance(self.ym1, self.yp1)
        return self.dy

